# 43. POL (スピンと偏極) の解説

最終更新: 2026-11-17

対象: `POL` / `RADPOL` フラグ、`Emittance[]` / `EMIT` の偏極の量 (偏極時間・スピンチューン・偏極ベクトル・スピン歳差行列・平衡偏極度)、`TrackParticles` のスピン追跡。
コードは `pysad/spin.py` (tradke のスピン部、spnorm、srequpol、spdepol、tchge の回転)、`pysad/emit6d.py`、`pysad/emit_extra.py`、`pysad/track6.py` (`_sprot`、`bsi`)、`pysad/ffscmd.py` (`fn_track`)。
試験は `pysad/tests/test_spin.py`、`test_spin_track.py`。参照はすべて SAD K64 の gs の出力。
短い要約は `docs/36_spin_csr.md` の §1 とスピン追跡の節にもある。K64 の癖は `docs/21_sad_bug_report.md` の B55〜B58。

## 1. 何ができるか

| やりたいこと | 方法 | 状態 |
|----|----|----|
| リングの偏極時間・スピンチューン・偏極ベクトル・スピン歳差行列 | `POL;` を付けて `Emittance[]` か `EMIT` | gs と一致 (§7) |
| 平衡偏極度 `Polarization` | 同上 | K64 が NaN を返す (B55)。pysad も NaN (§6) |
| 粒子ごとのスピン (歳差と Sokolov-Ternov の放射偏極) の追跡 | `POL; RAD;` を付けて `TrackParticles` (入力を 9〜11 列にする) | gs と 4e-12 以下で一致 (§4、§7) |
| 放射偏極だけ切る | `NORADPOL;` (`RADPOL` の既定は ON) | 追跡の Sokolov-Ternov 項だけに効く (§4) |

フラグの既定は `POL` が OFF、`RADPOL` が ON (`pysad/ffscmd.py` の `_FLAGS_DEFAULT`、K64 の `initb1.f` の `defflg`)。SADHelp の説明は「POL が ON のとき、RADPOL は追跡で Sokolov-Ternov 効果を入れる」。
POL が OFF のとき、偏極の量は K64 も 0 のまま (`PolarizationTime`・`SpinTune`・`Polarization` は 0、`PolarizationVector` は {0,0,0}、`SpinPrecessionMatrix` は単位行列)。`NominalSpinTune` は POL に関係なく `h0·g` (g = 0.001159652181280002 は電子の異常磁気能率、tfstk.f の gspin)。
`RADPOL` が `Emittance[]` の結果に効くかは、K64 の `temit.f` が `radpol0=radpol` を保存するだけで偏極の計算に使う箇所を見つけていない。pysad は `Emittance[]` では `RADPOL` を参照しない (追跡だけ)。

## 2. 使い方の例 (動く SAD スクリプト)

### 2.1 Emittance[] の偏極の量

`pysad/tests/test_spin.py` の平面リング (`PLANAR`、MOMENTUM 3 GeV、12 セル、`aγ` = 6.808) をそのまま使う。

```
MOMENTUM = 3GEV;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
CAVI RFC = (VOLT=1.5e6 HARM=64);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 12*CELL1 RFC);
FFS;
$FORM="18.14";
PageWidth = 1999;
USE RING;
CELL;
CALC;
RFSW; RAD; NOFLUC;
POL;
em = Emittance[];
Print["OUT ", {PolarizationTime, SpinTune, NominalSpinTune, Polarization, PolarizationVector}/.em];
EMIOUT; EMIT;
Exit;
```

gs の結果 (`test_planar_ring_matches_gs`): `PolarizationTime` = 123.8154010596169 s、`SpinTune` = 0.40407272052341、`NominalSpinTune` = 6.80814822478217、`PolarizationVector` = {0, 1, 0}、`Polarization` = NaN (§6)。
`EMIT` を `EMIOUT` の後に書くと、スピン歳差行列 (`sx sy sz` の 3×3) と、1 周の脱偏極ベクトル (`gx gy gz` の行、(x,px,y,py,z,pz) の 6 成分) が表に出る (書式は Fortran の `1p g15.7` まで一致、`test_emit_report_and_spin_precession_matrix_match_gs`)。
`SpinPrecessionMatrix/.em` は同じ 3×3 の行列。

`POL` を `NOPOL;` に替えると `{PolarizationTime, SpinTune}` は {0, 0} (`test_pol_off_leaves_spin_outputs_zero`)。

### 2.2 TrackParticles のスピン追跡

`pysad/tests/test_spin_track.py` の `pt1` (RAD ON・NOFLUC。平面リング + 空洞)。ビームは 9 列 `{x, px, y, py, z, dp, 角度, sy/|s|, flag}` で、10・11 列目は損失位置・損失周の予約 (§4)。

```
(平面リングは §2.1 と同じ。RING の定義まで)
FFS;
USE RING;
CELL;
CALC;
RFSW; RAD; NOFLUC;
POL;
em = Emittance[];
Print["EM ",{SpinTune,NominalSpinTune}/.em];
x0={1e-4,0};px0={0,0};y0={0,1e-4};py0={0,0};z0={0,0};dp0={0,1e-3};
r=TrackParticles[{1, {x0,px0,y0,py0,z0,dp0,{0,0},{0.1,0.2},{1,1},{1,1}}}, 20, 1, 3];
Print["R ",Length[r]," ",Dimensions[r[[2]]]];
Print["S ",r[[2]]];
```

戻りは `{位置, {x, px, y, py, z, dp, 角度, sy/|s|, flag}}`。2 粒子 (x と y を 1e-4 ずらした 1 番目と、dp=1e-3 の 2 番目)、スピンの初期値は (角度 0, sy/|s| = 0.1) と (0, 0.2)。3 周追跡した後の 7 列目 (`atan(sz, sx)`) と 8 列目 (`sy/|s|`) が gs と 4e-12 以下で一致する (`test_planar_ring_spin_tracking_matches_gs`)。
`RAD FLUC` (乱数あり) でも一致する (`pt3`、`pt13`)。`RAD NOFLUC` を `RAD FLUC` に替え、`SeedRandom` を gs と合わせれば同じ乱数の順になる (K64 の乱数の引き順まで写した)。

そのほかの動く例: `test_spin.py` の `COUPLED` (縦ずれ・回転した四極・六極・ROTATE 付き偏向のあるリング。偏極ベクトルが {-1.35e-4, 0.99954, 0.0305})、`test_mult_ring_matches_gs` (MULT)、`test_solenoid_ring_matches_gs` (SOL 領域)。
CSR などと一緒の記録は `docs/36_spin_csr.md`、作業の経緯は `docs/history/work_unimpl.md` (POL の節)。

## 3. 物理と K64 の式

記号は K64 (`photons.f` の `tradke`、module `tspin`) のまま。a = g = 1.15965e-3、γ = h1、`pr` = 1+δ。

### 3.1 1 ステップのスピン歳差 (tradke)

放射の 1 スライスごとに、運動量の変化 (px, py, pz の前後) からスライス内の磁場を作る。
- 磁場に比例する量 `xpx, xpy, xpz` (運動量の変化の外積) に、縦の磁場の補正 `bsir0` と SOL の `bzh` を足した `xpzb = xpz + (bsir0 + 2·bzh·al)·pr²`。
- スライス中央の運動方向 (pxh, pyh, pzh) への射影 `bp` から縦成分 bl = bp·(pxh, pyh, pzh)、横成分 bt = (xpx, xpy, xpzb)/pr² − bl。
- スピン方程式の回転ベクトル `g = (1 + aγ)·bt + (1 + a)·bl`。
- `srot` (3×9) の列 1〜3 をベクトル g のまわりに |g| 回す (Rodrigues の式、`xsincos` で 1−cos を精度良く)。列 4〜9 には座標微分 `dg` を `dg·transr` として積算する (微分は `spin.py` の `tradke_spin` で解析的に作る)。
- Sokolov-Ternov の傾き `gintd` に `(srot·ĝ)[:3]·sflc·anp·(bt·γ·p/al1)²` を足す (`sflc = 0.75·(ELRADI/FINEST)²`)。
- 最後に偏向角 `phi0` の分だけ `srot` を回す (`srot(1,:) = cp·srot(1,:) + sp·srot(3,:)`、`srot(3,:) = (-sp·srot(1,:) + srot(3,:))/cp`。2 つめは更新後の 1 行目を使う: B56)。

bt の z 成分の微分 `dbtz` は K64 のソースで `pr**2` ではなく `pr` で割っている。pysad はそのまま写している (`spin.py` の注記)。

### 3.2 座標回転と bsir0

- `tchge` / `tbrote` (ずれ・ロールのある偏向磁石の入口・出口): ROTATE に相当する回転で `srot(1:2,:)` を z 軸まわりに回し (`tchge_rotate_spin`)、ずれ・DROTATE・CHI2 は `tbrote` の微小回転 `rr` で更新する。`tbrote` の更新は列ごとに行を逐次更新する (rr·srot の行列積ではない: B56)。
- `bsir0` は縦の磁場の積分で、`tradke` の bl に入る。更新の場所は偏向磁石 (`tbendebody` の入口・出口の `bsi1·yi/rhob`、`bsi2` の項は捨てられる: B58)、`tmulke` (MULT の `Im(cx)/al`、最初のスライスで足し最後で引く。L=0 の薄い MULT では更新しない)、`tsetr0` (0 にリセット)。
- 偏向磁石の分割数は POL のとき `max(ndiv, 1 + int(max(|φ|, 1e-6)·h0·a/1e-2))` 以上になる (放射の結果も少し変わる。K64 も同じ。`dphipol = 1e-2`)。
- 傾いた四極・六極の入口・出口の回転 (`trot33`) は、pysad が実験室座標で `tradke` を計算するので入口と出口で打ち消し、書かなくても同じ結果になる (Emittance[] の側)。追跡側は §4。

### 3.3 偏極ベクトルとスピンチューン (spnorm)

`temit.f` の `tecalc` の最後 (`spin.py` の `spin_results`) で、1 周の行列 `srot(:,1:3)` を次の手順で処理する。
1. 列を直交化する (2 列目を規格化、1 列目から 2 列目の成分を引いて規格化、3 列目は外積)。
2. 回転の軸 (行列の反対称部分 `(R32−R23, R13−R31, R21−R12)`) を不変軸 n の初期値にする。大きさが 1e-4 未満なら実固有値の固有ベクトルを使う。
3. `(R − 1)·δ = n − R·n` の最小ノルム解 (`tsvd`、しきい値 1e-8) で軸を補正し、規格化、成分が負側に偏れば符号を反転する。これが `PolarizationVector`。
4. 軸に垂直な基底を作り、1 周での歳差角 `smu` を `atan2(-sm, cm)` で出す。`SpinTune` = smu/2π。
5. 減衰率 `sdamp = n·gintd`。`PolarizationTime = 1/sdamp/f_rev` (f_rev = ω0/2π は周回周波数。§5 の SOL の周長の癖が効く)。

### 3.4 平衡偏極度 (srequpol、spdepol)

`srequpol` は、規格化したモード座標でのビーム行列 (`beamn`)・減衰率・エミッタンス・チューンから 2 次までの摂動で脱偏極行列 `rm1` を作り (`spdepol`)、
`equpol = PST·sdamp / (sdamp − rm1[0,0])` (PST = 8√3/15 = 0.9238) を返す。`Polarization/.em` はその第 1 成分で、`Polarizations` は 3 成分 (残りは 0)。
SADHelp には `Polarization2`・`Polarization4`・`Polarization6` (2・4・6 次) もあるが、K64 のソースでは `spcalcres` を 2 次の 1 回だけ呼んでいる。pysad はそれに合わせ、`Polarization2` などの別の次数のキーは持たない (K64 のソースに `Polarization2` の文字列は見つからなかった。SADHelp にだけ記載がある。gs の `Emittance[]` も `Polarization2` を返さない (`Polarization2/.em` が未評価のまま)。pysad も同じ。2026-11-16 確認)。
`spdepol` の中で `exp(cb²−cφ)·cerfc(cb)` を作る。|cb| > 18 のとき `cerfc` が漸近式を使い、`exp(−cb²)→0` と `exp(cb²)→∞` の積が NaN になる (B55)。§6。

### 3.5 追跡のスピン (sprot)

`photons.f` の `sprot` は、粒子ごとに 1 ステップの磁場の向き (運動量の変化 ppx, ppy, ppz と `bsi`) から同じ g を作り、スピン (sx, sy, sz) を g のまわりに回す。
RADPOL が ON で、そのステップに放射エネルギー `ucde > 0` があれば、回転の前に Sokolov-Ternov の項 `dst = (PST − s·b̂t)·CFPD·(pp·p0/γ²)²·ucde/|bt|` (CFPD = 5√3/8) を bt の向きに足す。最後に偏向角の分だけ回す (sx の更新の後の sz は更新後の sx を使う: K64 のまま)。

## 4. TrackParticles のスピン追跡

- 入力と戻り: POL が ON のとき、ビームは 9〜11 列でなければならない (それ以外は `General::wrongtype`)。7 列目は歳差の角度、8 列目は `sy/|s|`、9 列目は flag、10・11 列目は損失位置・損失周。
  粒子のスピンは (sx, sy, sz) = (st·cos 角度, sy, st·sin 角度)、st = √(1−sy²)。戻りの 7 列目は `atan(sz, sx)`、8 列目は `sy/|s|` (tftrack.f の `tfsetparticles` / `ktfresetparticles`)。POL が OFF なら 7〜9 列のまま (`docs/08_tracking.md`)。
- 動く条件: **POL かつ RAD が ON のときだけ**。K64 は歳差を `tradk` の中で計算するので、RAD が OFF では高速追跡 (`TrackRun`) が走り、スピンは変わらないまま戻る (7・8 列は入力の値を正規化して返すだけ)。RAD が ON で POL なら `fn_track` は numba の高速追跡を使わず、`track6.py` の Python の追跡に切り替える (`ffscmd.py` が `NotImplementedError` で切り替え)。
  `RADPOL` が OFF なら Sokolov-Ternov の項だけが消え、歳差は動く。
- **損失の処理はない**: 高速追跡 (`TrackRun`) にしか損失の判定 (`LOSSAMPL`、`LOSSDZ`、flag) がないので、POL のスピン追跡では粒子は落ちない。flag 列は入力のまま返る。10・11 列目 (損失位置・周) は更新しない。
- `track6.py` の構成: `_Particle.s` (スピン)、`_Particle.bsi` (縦の磁場の積分)、`_sprot` (§3.5)、`_tradk` の中で呼ぶ (FLUC で光子が 0 個の場合も歳差だけ呼ぶ: K64 の `tradkfn` の elseif)、`_spin_rot_z` (四極の ROTATE の chi3 回転: `tsolrot` の calpol)。
- bsi の入る場所 (K64 の対応する関数): 四極 (`tsolqur` の各スライスの両端)、BEND (`tbendcore` の入口・出口の y/rhob、K1 のある BEND は `tbendi` の各スライスの `akxn·xi·yi`)、MULT (`tmulti` の `Im(.5·cx·cx1²)/al`)、`tdrift` の ak0 付き、SOL 領域内のドリフトと四極 (`bz·al`)。

## 5. 要素ごとの扱いと K64 の癖

### 5.1 SOL 領域

- `Emittance[]`: SOL (BOUND≠0 の領域) を含むリングで偏極時間・スピンチューン・偏極ベクトルが gs と一致する (`emit_extra._tturne1` が `_solmap` の領域を通す。`test_solenoid_ring_matches_gs`、領域内に四極 QT を含む)。結果 (偏極時間・スピンチューン・偏極ベクトル) が一致したことだけが確認済みで、領域内の bsi の入れ方を要素ごとに分けた照合はしていない。
- **K64 は SOL (BOUND≠0) の長さを周長に入れない** (`emit6d._emit_circ`。PLANAR リングに L=0.3・L=0.6・BZ=0 の SOL を足しても、周長 45.6 のときと偏極時間が 15 桁一致)。周回周波数・`MomentumCompaction`・偏極時間に効く。pysad も同じにしてあり、追跡の `circ` も `_emit_circ` を使う。
- SOL に GEO がないと、K64 の S 座標が途切れて周長がさらに狂う (偏極時間が約半分)。pysad は再現しない (未対応)。
- TrackParticles: SOL 領域内のドリフトは gs と一致 (`pt8`、BZ=0.05、MOMENTUM 3 GeV)。
- **領域内の MULT・BEND・CAVI は未対応** (Emittance・追跡とも)。**領域内の四極のスピン追跡は gs の参照が食い違い、未照合** (§8)。

### 5.2 MULT

K1・K2・K3 とずれ (DX・DY) のある MULT のリングは Emittance[] で gs と 1e-11 で一致 (`test_mult_ring_matches_gs`。K3=300 まで)。`tmulke` の `bsir0` の更新 (`kick_cx`、`_tmulte`) を写した。L=0 の薄い MULT では更新しない (以前は 0 割りになっていた: `test_thin_mult_ring_spin_tune_matches_gs`)。
追跡も K1 K2 K3 の MULT と薄い MULT で 4e-12 以下 (`pt5`、`pt7`、`pt7b`)。ANGLE≠0 の MULT (`tmultae`) の POL スピンも gs と一致 (偏極時間 14 桁、スピンチューン 1e-14。偏極ベクトルは平面リングで 0,1,0。2026-11-16)。

### 5.3 空洞

偏向のない空洞は、Emittance 側で放射と偏向がないので `srot` は変わらないはず、というのがコード上の理解。`emit6d` の空洞 (`_tcave`) には `tchge_spin` の呼び出しがなく (BEND の `_tbende` だけにある)、DX・DY・ROTATE のある空洞の入口・出口の回転は `srot` に反映しない。入口と出口の回転が打ち消すので結果は同じとみられるが、ROTATE+DX+DY 付き空洞のスピンは結合リングで gs と一致 (偏極時間 13 桁、スピンチューン 1e-14、偏極ベクトル 1e-13。2026-11-16。自動試験には入れていない)。
RF 空洞を含むリングそのもの (`RFC`) は、すべての試験の例に入っていて一致している。

### 5.4 K64 の癖 (docs/21)

| 番号 | 内容 | pysad |
|----|----|----|
| B55 | `Polarization` が NaN (`spdepol` の `cerfc`) | 同じ式で NaN (§6) |
| B56 | `tbrote`、`serot`、`tradke` 末尾の `srot` 更新が `do concurrent` の中で行を逐次更新 (2 行目は更新後の 1 行目を使う) | そのまま写す (`serot`・`tradke` 末尾・`tbrote` とも gs で確認。`tbrote` は BEND の DROTATE+DX で、偏極時間 14 桁、スピンチューン 1e-14、偏極ベクトル 1e-14 で一致) |
| B57 | `tadddsg` の呼び出しが dgz に dgx を渡す (結果に影響しない) | `dsg` を持たない |
| B58 | `bsi2` の項が出口で捨てられる (縦のオフセットのある BEND) | 捨てる |
| SOL | SOL の長さを周長に入れない | 同じ (§5.1) |

そのほか: `dbtz` の割り算 (`pr` で割る)、`sprot` の最後の `sz` 更新が更新後の `sx` を使う (K64 のまま)。

## 6. Polarization (平衡偏極度) が NaN になる件

B55 を参照。gs で再現し、pysad も同じ式で NaN になる。試した 8 つのリング (平面でも非平面でも) すべてで NaN。
`cb = −(e^{−cφ}−1)/a` (`a = √2·ssprd`、ssprd はスピン軌道結合 ~0.05 など) が通常 20 以上になり、漸近式の `exp(−cb²)` と `exp(cb²)` の積が 0×∞ になるため。
参考: `erfc·exp` を桁あふれしない形 (erfcx) に書き直すと有限になるが、値は 4e-7〜8e-7 (つまり 4e-5〜8e-5 %) で、92% 級の理論値からは遠く、`spdepol` の式か単位に誤りがある疑いがある (推定。K64 の式の正しさは未確認)。
**未確認**: gs で `Polarization` が有限になる例は不安定な光学の 1 つだけで、pysad はその光学の閉軌道を作れない (正規形が不安定) ので、`cerfc` の |cb| ≤ 18 側は gs で照合できていない。
実用上は `PolarizationTime` (= 1/ΣT の式で、有限で正しい) を使う。

## 7. gs との一致の度合い

試した格子: エネルギー 3.0〜3.08 GeV (aγ = 6.8〜6.99)、縦ずれ DY 0〜2 mm、ROTATE 0〜0.03 rad。

| 量 | 平面リング | 縦ずれ・傾き・六極のあるリング |
|----|----|----|
| PolarizationTime | 全桁 (123.8154010596169 s) | 13 桁 (130.35269523782 s) |
| SpinTune | 13 桁 | 12 桁 |
| PolarizationVector | 12 桁 | 10 桁以上 |
| SpinPrecessionMatrix・脱偏極ベクトルの表 | 一致 (書式まで) | 一致 |
| Emittances・DampingRate (POL 付き) | 11〜12 桁 | 11〜12 桁 |
| TrackParticles の戻り (7・8 列) | 4e-12 以下 | 4e-12 以下 (NOFLUC・FLUC) |

- MULT のリング、SOL (ドリフトと四極のある領域) のリングは偏極時間・スピンチューン 1e-11、偏極ベクトル 1e-10 (`atol`) で一致。追跡は平面リング、ずれ・回転のあるリング、MULT (K1 K2 K3)、K1 のある BEND、SEXT、薄い MULT、SOL 領域内のドリフトで 4e-12 以下 (試験の許容は `tol=1e-8` 相対)。
- すべて MOMENTUM 3 GeV 前後、周長は十数 m〜数十 m の小さな試験リング (周長に SOL の長さは入らない)。実機の格子での照合はしていない。
- 試験は 13 件が通る (`pytest pysad/tests/test_spin.py pysad/tests/test_spin_track.py`、2026-11-16 に確認、約 12 秒)。

## 8. 限界・未確認

- **SOL 領域内の MULT・BEND・CAVI は未対応**。BEND は K64 の `tsol.f` (`tdrift` を呼ぶ) どおりに直訳したが、gs が偏向をほぼ無視して (BZ=0.5 T の領域で ANGLE=0.1 でも位置・角度が 1e-5 しか変わらない) 合わず取り下げた。原因は未確認 (`docs/history/work_unimpl.md`)。
- **SOL (BZ≠0) 領域内の四極のスピン追跡は gs の参照が食い違い、未照合**。SOL + 四極 + RAD の追跡そのものが、保存済みの参照 (`bench_level9.sad`、2026-10-14) と現在の gs (2026-10-19 ビルド) で別の値になる (x が −0.0027867898 と −0.0014539780)。NORAD は両方で一致する。どちらが正しいかは断定できないが、2026-11-16 の切り分けでは現在の gs 側が SOL 領域内の四極 + RAD を誤って扱う疑いが強い (`docs/44_sol_guide.md` §7 の限界 1)。
- BEND の `F1=0.1 FRINGE=3` を RAD ON で追跡すると gs と合わない (相対 0.7)。L≠0 の CAVI に DX・DY を付けると 1e-3 ずれる。どちらもスピンなしでも出る差で、未調査。スピン追跡の一致の範囲はこの 2 つを含まない。
- K3 が強くずれの大きい MULT リング (K3=300、ずれ 1e-4) は、pysad の `CALC` の閉軌道が gs と別の解 (y = −0.033) になり、`Emittance[]` が NaN になる。gs は小さい軌道で通る。`qcod` の初期値・後退の違いの可能性があるが未解決。
- `Polarization` は NaN (§6)。`Polarization2/4/6` は実装していない。
- ROTATE 付き空洞の `Emittance[]` のスピン回転は照合していない (§5.3)。
- 追跡は損失の処理なし、numba の高速追跡 (`TrackRun`) を使わないので遅い (Python の追跡)。
- `tbrote` の更新 (ずれ・DROTATE・CHI2 のある偏向磁石) は gs で確認していない (B56)。
- `Emittance[]` が SOL を含むときの偏極の量は上記の SOL の例 (四極 1 つ + ドリフト) だけ照合している。SOL の中の BEND・MULT・CAVI がないリングに限る。
- スロー経路 (`PYSAD_NO_ELKERN=1` 相当) の `_tquade` が、ずれだけで回転のない四極のヤコビアンをずらす前の座標で評価していた誤りは直した (POL は numba の速い経路を使えないので εy に 7e-5 出ていた)。

## 9. 関連文書

- `docs/36_spin_csr.md`: POL と CSR の要約 (§1 とスピン追跡の節)
- `docs/29_sol_optics.md`: SOL の Twiss・Emittance・追跡
- `docs/07_emittance.md`、`docs/08_tracking.md`: Emittance と TrackParticles の一般説明 (POL の記述は更新が必要な箇所がある)
- `docs/21_sad_bug_report.md`: K64 の癖 B55〜B58
- `docs/history/work_unimpl.md`: POL・SOL・追跡の作業記録と未解決の差
- `docs/45_radlight_guide.md`: RADLIGHT (`TrackParticles` が軌道を返す)
