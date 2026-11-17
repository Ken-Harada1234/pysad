# 36. POL (スピン) と CSR

最終更新: 2026-11-17 (現状のスナップショット。経緯は `history/` にある)

`POL;` を付けた `Emittance[]` / `EMIT` のスピン (偏極時間・スピンチューン・偏極ベクトル・スピン歳差行列) と、`TrackParticles` のスピン追跡 (§3)、CSR のウェイク追跡の関数 (§2)、RADLIGHT (§4)。
K64 の `photons.f` (tradke、tspin モジュール)、`tchge.f`、`temit.f`、`csrtrack.f`、`csroy.f`、`csrtest.f` の直訳。
照合は SAD K64 (gs) の出力。コードは `pysad/spin.py`、`pysad/csr.py`、テストは `pysad/tests/test_spin.py` (8 本)、`test_spin_track.py` (5 本)、`test_csr.py` (7 本。うち 1 本は `PYSAD_SLOW=1`)、`test_radlight.py` (1 本)。

## 1. POL

### 1.1 使い方

```
FFS; ... CELL; CALC; RFSW; RAD; NOFLUC;
POL;
em = Emittance[];
PolarizationTime/.em   SpinTune/.em   PolarizationVector/.em   SpinPrecessionMatrix/.em   Polarization/.em
EMIOUT; EMIT;          (* スピン歳差行列と 1 周の脱偏極ベクトルも表で出る *)
```

`POL` を付けない場合は従来どおり 0 (K64 も params=0 のまま)。

### 1.2 何を計算するか

- `srot(3,9)`: 列 1〜3 がスピン歳差の 1 周行列 (初期値は単位行列)、列 4〜9 が (x,px,y,py,z,pz) の初期座標に対する 1 周のスピン微分 (脱偏極ベクトル)。
- 放射の 1 ステップ (`tradke`) ごとに、スライス内の磁場の横成分 (bt)・縦成分 (bl)・スピン方程式の回転ベクトル g=(1+aγ)bt+(1+a)bl とその座標微分 dg を作り、`srot` を回して微分を積算する (spin.py の `tradke_spin`)。
  放射の傾き `gintd` (Sokolov-Ternov の向き) も積算し、1/(n·gintd)/f0 が偏極時間 (秒)。
- `tchge` の座標回転 (ROTATE の付いた偏向磁石) で `srot` を回す。`tbrote` の微小回転も写した (gs 未照合)。
- 偏向磁石の分割数は POL のとき `1+int(|φ|·h0·a/0.01)` 以上になる (放射の結果も少し変わる。K64 も同じ)。
- 傾いた四極・六極 (tsolrote 側の `trot33`) は、pysad が実験室座標で `tradke` を計算するので回転しなくても同じ結果になる (回転は入口と出口で打ち消す)。
- `spnorm`: `srot` を直交化して不変軸 n (PolarizationVector) とスピンチューン、減衰率を求める。`srequpol` → `spdepol`: 2 次の摂動で脱偏極行列を作り、平衡偏極度 equpol を出す。

### 1.3 gs との一致

| 量 | 平面リング | 縦ずれ・傾き・六極のあるリング |
|----|----|----|
| PolarizationTime | 全桁 (123.8154010596169 s) | 13 桁 (130.35269523782 s) |
| SpinTune | 13 桁 | 12 桁 |
| PolarizationVector | 12 桁 | 10 桁以上 |
| SpinPrecessionMatrix・脱偏極ベクトルの表 | 一致 (書式まで: Fortran の 1p g15.7) | 一致 |
| Emittances・DampingRate (POL 付き) | 11〜12 桁 | 11〜12 桁 |

エネルギー 3.0〜3.08 GeV (aγ=6.8〜6.99)、縦ずれ DY 0〜2 mm、ROTATE 0〜0.03 rad の格子を試した。

### 1.4 限界・未確認

- `Polarization` (平衡偏極度) は K64 が NaN (B55) で、pysad も同じ式で NaN になる。有限値になる gs の例は不安定な光学の 1 つだけで、pysad はその光学の閉軌道を作れない (正規形が不安定で失敗) ため、`spdepol` の有限値側 (|cb|≤18 の cerfc) は gs で確認できていない。
  erfcx で桁あふれを避けると 4e-7〜8e-7 の有限値になる (この値も理論値から遠く、K64 の式か単位の疑い)。
- Emittance[] / EMIT のスピンは次まで対応し gs と一致する (`test_spin.py`): MULT (K1・K2・K3・ずれ。`tmulke` の bsir0 更新も入れた。L=0 の薄い MULT を含む)、空洞 (試験のリングはどれも `CAVI RFC=(VOLT=1.5e6 HARM=64)` を含む)、
  SOL (`emit_extra._tturne1` が SOL 領域 `_solmap` を通す。GEO は 1 つだけ、`test_solenoid_ring_matches_gs`)。K64 は SOL (BOUND≠0) の長さを周長に入れず (omega0・MomentumCompaction・偏極時間に効く)、`emit6d._emit_circ` で同じにした
  (PLANAR + SOL の偏極時間が、SOL を除いた周長 45.6 のときと 15 桁一致)。GEO のない SOL では K64 の S 座標が途切れて周長がさらに狂う (偏極時間が約半分) が、pysad は再現しない。
- 2026-11-16 夜に追加で gs と一致を確認 (CALC・POL の Emittance): ROTATE+DX+DY 付きの空洞 (結合リング。偏極時間 13 桁、スピンチューン 1e-14、偏極ベクトル 1e-13)、ANGLE≠0 の MULT (`tmultae`、docs/28。偏極時間 14 桁、スピンチューン 1e-14。平面リングなので偏極ベクトルは 0,1,0)、BEND の DROTATE+DX (`tbrote`。偏極時間 14 桁、スピンチューン 1e-14、偏極ベクトル 1e-14)。いずれも 1 つのリングでの照合で、広い格子ではない。
- 未対応・未確認: `MAP` のスピン、SOL の領域内の MULT・BEND・CAVI (領域内の要素そのものが未対応: docs/29)。
  K3 が強くずれの大きい MULT のリング (K3=300、ずれ 1e-4) は、pysad の `CALC` (qcod) の閉軌道が gs と別の大きい値 (y=-0.033) になり、Emittance が NaN になる (gs は小さい軌道で通る)。qcod の初期値・後退の違いの可能性で未調査。
  縦のオフセットのある偏向磁石の bsi の項 (`bsi2`) は K64 でも捨てられる (B58)。
- スロー経路 (`PYSAD_NO_ELKERN=1` 相当) の `_tquade` は、ずれ (DX/DY) だけで回転のない四極のヤコビアンを、ずらした後の座標で評価する (ずらす前の座標で評価すると、POL は numba の速い経路を使えないので εy に 7e-5 の誤差が出る)。
- K64 の癖は B55〜B58 (docs/21)。

## 2. CSR

`CSRTrack`、`CSRHaissin` (csrtrack.f)、`CSRSetupOY`、`CSROYMatrix` (csroy.f)、`CSRMatrix`、`CSRConvert` (csrtest.f の一部)。実FFT (`tcftr`・`trftr`・`tftrr`) と自然スプライン (`spline1`) も写した。

| 関数 | 状態 |
|----|----|
| `CSRTrack[z, dp, zl, {dz, ndiv, inten, r56, r65, frf, phirf, damp, sigmae}, nturn]` | gs と 11 桁一致 (乱数は ffs と共通の SadRandom、z・dp はその場で更新、発散後は 0、RF あり・なし) |
| `CSRHaissin[zl, {dz, inten, r65, sigz}]` | gs と 11 桁一致 (Newton 反復。tsolvg は pysad/tsvd.py) |
| `CSRSetupOY[rho, zl, {dz, nj, nphi, inten, r56, r65, zspan}]` | gs と J・ω 1e-9、rhoj 1e-9 (周回の離散計算のため)。dpjn の記憶域の重なり (B59) を再現 |
| `CSROYMatrix[...]` | ソースどおりに翻訳 (gs が wrongtype を返すので未照合: B60) |
| `CSRMatrix[pipe, {dx,dy,f0..f3,tilt}, omega, rho, center, mode]` | mode 22・33・23・32・110・220・330 は gs と全要素一致。mode 1000 (らせん断面 `csrhelmat`) は実装済みで 2 断面が相対 1e-13 で一致 (§2.1) |
| `CSRConvert[v or m, pipe]` | gs と一致 (K64 の tfcsrvec と同じ) |
| `CSRInit[{{dx,dy}, center, rho, pipe}]` | 実装済み、gs と相対 1e-10 (§2.1) |

### 2.1 CSRInit・mode 1000・CSR.n の流れ (実装済みの範囲と未照合の範囲)

- `CSRInit[{{dx,dy}, center, rho, pipe}]` は実装済み (csrtest.f の csrphimesh・csr2dpoisson。境界の電荷は `tsolva` = `tsvd(.false.)` の最小ノルム解で、pysad の `tsvd` を使う。
  `csr2dpoisson` の電場は 3 回の繰り返し)。gs と全要素の和・最大値・個別の値が相対 1e-10 で一致 (pipe 8 列と 5 列の 2 例)。実数 1 つの dx は K64 のバグ
  (dx を読まず 0 のまま) で常に wrongtype。試験は `test_csr.py` の `test_csrinit_matches_gs`。
- `CSRMatrix` mode 1000 (らせん断面。`csrhelmat` + `csr2dpoisson1`) は実装済み。境界の点列を中点で 2 倍にして境界要素法で各セルの電場への応答 (2nm×2nm) を求める。
  `tsolvm(svd=.true.)` は最小ノルム解 (特異値 < 最大×1e-4 を捨てる) として numpy の SVD で解いた。K64 は `phix(0:nphimax)` を `csrphimesh` の `phix(*)` に渡すため、
  点の添字が 1 つずれる (これも写した)。gs と {2,2,2} と {5,5,5,5,5,5} の 2 つの断面で全要素の和・個別の値が相対 1e-13 で一致。{4,4,4,4} などは gs が SIGABRT/SIGSEGV で落ちる (参照なし)。
- `CSR.n` の流れが動く: `CSR[…]` の生成、`MakePipe[]`、`CSRData` (CSRMatrix → RealEig → CSRInit → CSRConvert)、`Ez`、`MakeZL` (周波数の表からインピーダンス)。
  8 列のメッシュの `CSRData` は固有値が gs と相対 1e-15、`Ez` が相対 5e-10 で一致 (固有ベクトルは縮退で並びが違うが物理量は同じ)。
  `CSR[NParallel->1, Nk->4, BunchLength->3e-4, BendingRadius->10, PipeWidth->0.04, Center->0.02]@MakeZL[{{10.,1.}}, 2000., 6000., 3]` の ZL の表 (3 周波数) は
  gs と相対 1e-10 以内で一致 (pysad で約 7 分、試験は `PYSAD_SLOW=1` で動く `test_csr_makezl_matches_gs`)。
  動かすために直したもの: `Do[…, {Infinity}]`・`{i, Infinity}` (Break するまで続ける)、`OpenRead["!コマンド"]` (コマンドの出力を読む。`top` のない環境 (Windows) では
  `FreeMem` 用に空き記憶領域を旧形式の top の行で作る)、開けないファイルは `$Failed`。
- gs の参照を取るには、WSL に `top` を入れるか、PATH に `top` の代用 (`Mem: …k total, …k used, …k free, …k buffers` を出すスクリプト) を置き、`OSTYPE=linux` を設定する
  (非対話の WSL では OSTYPE が空で、`FreeMem` の `Switch` が未評価になり OpenRead が失敗する)。
- `Packages/CSR.n` (1088 行) は `pkgload.py` の許可リストにあり、読み込める (描画の `ListPlot.n`・`Plots.n` なども許可リストにある)。
- **未照合 (gs と比べていない)**: `CSROYMatrix` (gs が wrongtype を返す: B60)、`CSR.n` の `ScanTrack`・`OYStability` などの長い流れ、実用の大きさのメッシュの `MakeZL` (試験の 3 周波数・2000〜6000 の範囲以外)、
  mode 1000 の {4,4,4,4} など gs が落ちる断面 (参照なし)、実数 1 つの `CSRInit` (K64 のバグで常に wrongtype)。`CSR.n` の `SkinDepth`・`PathDifference`・`StandardKList`・`VariableMesh`・`ZThre` は gs と一致 (docs/history/work_unimpl.md の記録。試験は MakeZL だけ)。

- 2026-11-16 夜にソースで確認: `CSR.n` の `TrackAnimate` は `ListPlot` を呼んだあと `System["import -silent -window … -format -jpeg …"]` (ImageMagick) と `System["mkdir …"]` を呼ぶので、Windows の pysad (System は実行しない) では画像は作れない。
  `Wait[mem]` は `FreeMem[]` (top) と `Sleep` を使う。gs はこの環境で top が無く FreeMem が OpenRead エラーになる。pysad の `FreeMem` は Python (psutil) で代替して 1.0 を返す。並列を増やすと速くなるかは確かめていない (pysad は fork しない)。

K64 の癖は B59・B60。

## 3. スピン追跡 (TrackParticles、POL)

POL と RAD が ON のとき、`TrackParticles` は粒子ごとのスピンを追う (K64 では tradk の中の `sprot` で歳差と放射偏極を計算するので、RAD OFF ではスピンは変わらない)。
実装は `pysad/track6.py` (Python の追跡。高速追跡 TrackRun は使わないので、粒子の損失 (flag) の処理はない)。

- 入力は 9〜11 列 `{x,px,y,py,z,dp, 角度, sy/|s|, flag}`、戻りの 7・8 列は `atan(sz, sx)` と `sy/|s|` (tftrack.f の `tfsetparticles` / `ktfresetparticles`)。スピン (sx, sy, sz) = (st·cos 角度, sy, st·sin 角度)、st = √(1−sy²)。
- bsi (磁場の縦成分の積分) が入る場所: 四極 (tsolqur の各スライスの両端)、BEND (tbendcore の入口・出口、K1 のある BEND は tbendi の各スライス)、MULT (tmulti)、tdrift の ak0 付き、SOL 内のドリフトと四極 (bz·al)。四極の ROTATE は tsolrot の chi3 のスピン回転。
- gs と一致 (4e-12 以下、NOFLUC・FLUC): 平面リング、ずれ・回転のあるリング、MULT (K1 K2 K3)、K1 の BEND、SEXT、薄い MULT、SOL 領域内のドリフト。`test_spin_track.py`。
- 限界: SOL (bz≠0) 内の四極のスピンは gs と照合できていない (SOL + 四極 + RAD の追跡そのものが、現在の gs (2026-10-19 ビルド) と保存済みの参照 (bench_level9、10-14) で食い違う。11-16 の単一粒子の調べでは、gs が RAD ON・BZ≠0 の領域内で四極の力を入れていない疑いが強く、pysad は物理と保存済みの参照に合う。ソースのどこかは未特定: docs/29、docs/history/work_unimpl.md)。
  BEND の `F1=0.1 FRINGE=3` (RAD ON、相対 0.7) と L≠0 の CAVI の DX・DY (相対 1e-3) は、スピンなしでも追跡が gs と合わない (未調査)。DX・DY・ROTATE だけの BEND は 1e-13 で一致する。
- 損失 (flag) の処理はない (高速追跡 TrackRun を使わないため、粒子は落ちない)。スピンは POL かつ RAD が ON のときだけ動く。K64 の `tradk` の中の計算なので、RAD OFF ではスピンは変わらない。ほかに、SOL 領域内の MULT・BEND・CAVI は追跡そのものが未対応 (docs/29)。

## 4. RADLIGHT

`TrackParticles` は RADLIGHT のとき RAD・TRPT を ON にし、`{{ls または ld, 座標}, 軌道}` を返す。軌道は粒子ごとに {t, x, y, z} の表 (`RadiationField` の入力)。
K64 の tlinit / tlstore / tlresult の直訳 (`ffscmd._RadLight`、gs と差 0: `test_radlight.py`)。K64 の癖: 記録するのは区間 (tturn) の始めと、L≠0 の BEND 以外の要素の後だけで、BEND の中は記録しない。
横方向の座標は軌道に入らない (tlinit が geo0 の軸 1・2 を設定しない)。位置は straight (dphi=0) の積み上げで BEND で向きが変わらない。
`RadiationField`・`RadiationSpectrum` は実装済み。BEND の中を記録しない K64 の仕様のため、SAD 付属の例の流れ (B を ANGLE 0 の STEER にする) でないと意味のある場にならない。
2026-11-16 夜の確認:
- 軌道の「進行方向」は局所 z 軸で、`geo[0]=[[0,0,1,0],[-1,0,0,0],[0,-1,0,0]]` から局所 z 軸の大域成分は (1,0,0) なので、進行方向の距離は大域 x (軌道の 2 行目) に出る。
- RADLIGHT → RadiationField → FindPeak → RadiationSpectrum の流れ (B を ANGLE=0 K0=0.1 の STEER にして 1 周): pysad は最後まで通る (軌道 8 点 {4,8}、f {10,8}、ピーク位置 3、スペクトル {3,1001})。gs は RadiationField で SIGSEGV (B49) のため通しの比較はできない。軌道 (TrackParticles の戻り) の部分だけ gs と一致。
- 限界: 複数区間 (複数周) の軌道は区間の始めの点が重なる (t が同じ点が 2 つ) ため、pysad の `RadiationField` は `x must be strictly increasing sequence` で失敗する。gs も RadiationField が SIGSEGV で落ちるので比較できない。
未確認: 大きなリング・多数粒子での規模。解説は docs/45_radlight_guide.md。

解説: docs/43_pol_guide.md (POL・スピン)、docs/40_csr_guide.md (CSR)、docs/45_radlight_guide.md (RADLIGHT)
