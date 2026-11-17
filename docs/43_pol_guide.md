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
| 平衡偏極度 `Polarization` | 同上 | K64 が NaN を返す (B55)。**k64 モードは NaN、v2 モードは有限値** (§6、§8a) |
| 粒子ごとのスピン (歳差と Sokolov-Ternov の放射偏極) の追跡 | `POL; RAD;` を付けて `TrackParticles` (入力を 9〜11 列にする) | k64 モード: gs と 4e-12 以下で一致 (§4、§7)。v2 モード: 一部の歳差の式を直す・領域内や損失も扱う (§8a) |
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
- **損失の処理 (k64 モード)**: 高速追跡 (`TrackRun`) にしか損失の判定 (`LOSSAMPL`、`LOSSDZ`、flag) がないので、k64 モードのスピン追跡では粒子は落ちない。flag 列は入力のまま返る。10・11 列目 (損失位置・周) は更新しない。**v2 モードは track6 の経路でも落とす** (§8a)。
- `track6.py` の構成: `_Particle.s` (スピン)、`_Particle.bsi` (縦の磁場の積分)、`_sprot` (§3.5)、`_tradk` の中で呼ぶ (FLUC で光子が 0 個の場合も歳差だけ呼ぶ: K64 の `tradkfn` の elseif)、`_spin_rot_z` (四極の ROTATE の chi3 回転: `tsolrot` の calpol)。
- bsi の入る場所 (K64 の対応する関数): 四極 (`tsolqur` の各スライスの両端)、BEND (`tbendcore` の入口・出口の y/rhob、K1 のある BEND は `tbendi` の各スライスの `akxn·xi·yi`)、MULT (`tmulti` の `Im(.5·cx·cx1²)/al`)、`tdrift` の ak0 付き、SOL 領域内のドリフトと四極 (`bz·al`)。

## 5. 要素ごとの扱いと K64 の癖

### 5.1 SOL 領域

- `Emittance[]`: SOL (BOUND≠0 の領域) を含むリングで偏極時間・スピンチューン・偏極ベクトルが gs と一致する (`emit_extra._tturne1` が `_solmap` の領域を通す。`test_solenoid_ring_matches_gs`、領域内に四極 QT を含む)。結果 (偏極時間・スピンチューン・偏極ベクトル) が一致したことだけが確認済みで、領域内の bsi の入れ方を要素ごとに分けた照合はしていない。
- **K64 は SOL (BOUND≠0) の長さを周長に入れない** (`emit6d._emit_circ`。PLANAR リングに L=0.3・L=0.6・BZ=0 の SOL を足しても、周長 45.6 のときと偏極時間が 15 桁一致)。周回周波数・`MomentumCompaction`・偏極時間に効く。pysad も同じにしてあり、追跡の `circ` も `_emit_circ` を使う。
- SOL に GEO がないと、K64 の S 座標が途切れて周長がさらに狂う (偏極時間が約半分)。pysad は再現しない (未対応)。
- TrackParticles: SOL 領域内のドリフトは gs と一致 (`pt8`、BZ=0.05、MOMENTUM 3 GeV)。
- **領域内の MULT・BEND・CAVI は k64 モードでは未対応** (Emittance・追跡とも)。v2 モードは追跡のスピンが通る (§8a。Emittance[] の POL の量は領域内だけ K64 の式のまま)。**領域内の四極のスピン追跡は gs の参照が食い違い、未照合** (§8)。

### 5.2 MULT

K1・K2・K3 とずれ (DX・DY) のある MULT のリングは Emittance[] で gs と 1e-11 で一致 (`test_mult_ring_matches_gs`。K3=300 まで)。`tmulke` の `bsir0` の更新 (`kick_cx`、`_tmulte`) を写した。L=0 の薄い MULT では更新しない (0 割りになるため: `test_thin_mult_ring_spin_tune_matches_gs`)。
追跡も K1 K2 K3 の MULT と薄い MULT で 4e-12 以下 (`pt5`、`pt7`、`pt7b`)。ANGLE≠0 の MULT (`tmultae`) の POL スピンも gs と一致 (偏極時間 14 桁、スピンチューン 1e-14。偏極ベクトルは平面リングで 0,1,0。2026-11-16)。

### 5.3 空洞

偏向のない空洞は、Emittance 側で放射と偏向がないので `srot` は変わらないはず、というのがコード上の理解。`emit6d` の空洞 (`_tcave`) には `tchge_spin` の呼び出しがなく (BEND の `_tbende` だけにある)、DX・DY・ROTATE のある空洞の入口・出口の回転は `srot` に反映しない。入口と出口の回転が打ち消すので結果は同じとみられるが、ROTATE+DX+DY 付き空洞のスピンは結合リングで gs と一致 (偏極時間 13 桁、スピンチューン 1e-14、偏極ベクトル 1e-13。2026-11-16。自動試験には入れていない)。
RF 空洞を含むリングそのもの (`RFC`) は、すべての試験の例に入っていて一致している。

### 5.4 K64 の癖 (docs/21)

| 番号 | 内容 | pysad |
|----|----|----|
| B55 | `Polarization` が NaN (`spdepol` の `cerfc`) | k64: 同じ式で NaN (§6)。v2: 有限値 (§8a) |
| B56 | `tbrote`、`serot`、`tradke` 末尾の `srot` 更新が `do concurrent` の中で行を逐次更新 (2 行目は更新後の 1 行目を使う) | k64: そのまま写す (`serot`・`tradke` 末尾・`tbrote` とも gs で確認。`tbrote` は BEND の DROTATE+DX で、偏極時間 14 桁、スピンチューン 1e-14、偏極ベクトル 1e-14 で一致)。v2: `tbrote` の更新を (I+rr)·srot の行列積にする (§8a) |
| B57 | `tadddsg` の呼び出しが dgz に dgx を渡す (結果に影響しない) | `dsg` を持たない |
| B58 | `bsi2` の項が出口で捨てられる (縦のオフセットのある BEND) | k64: 捨てる。v2: 入口・出口の縦磁場の積分を入れる (§8a) |
| SOL | SOL の長さを周長に入れない | 同じ (§5.1) |

そのほか: `dbtz` の割り算 (`pr` で割る)、`sprot` の最後の `sz` 更新が更新後の `sx` を使う (K64 のまま)。

## 6. Polarization (平衡偏極度) が NaN になる件 (k64 モード。v2 モードは §8a で有限値)

B55 を参照。gs で再現し、**k64 モードの** pysad も同じ式で NaN になる。試した 8 つのリング (平面でも非平面でも) すべてで NaN。
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

- **SOL 領域内の MULT・BEND・CAVI は k64 モードでは未対応** (v2 モードは §8a)。BEND は K64 の `tsol.f` (`tdrift` を呼ぶ) どおりに直訳したが、gs が偏向をほぼ無視して (BZ=0.5 T の領域で ANGLE=0.1 でも位置・角度が 1e-5 しか変わらない) 合わず取り下げた。原因は未確認 (`docs/history/work_unimpl.md`)。
- **SOL (BZ≠0) 領域内の四極のスピン追跡は gs の参照が食い違い、未照合**。SOL + 四極 + RAD の追跡そのものが、保存済みの参照 (`bench_level9.sad`、2026-10-14) と現在の gs (2026-10-19 ビルド) で別の値になる (x が −0.0027867898 と −0.0014539780)。NORAD は両方で一致する。どちらが正しいかは断定できないが、2026-11-16 の切り分けでは現在の gs 側が SOL 領域内の四極 + RAD を誤って扱う疑いが強い (`docs/44_sol_guide.md` §7 の限界 1)。
- BEND の `F1=0.1 FRINGE=3` を RAD ON で追跡すると gs と合わない (相対 0.7)。L≠0 の CAVI に DX・DY を付けると 1e-3 ずれる。どちらもスピンなしでも出る差で、未調査。スピン追跡の一致の範囲はこの 2 つを含まない。
- K3 が強くずれの大きい MULT リング (K3=300、ずれ 1e-4) は、pysad の `CALC` の閉軌道が gs と別の解 (y = −0.033) になり、`Emittance[]` が NaN になる。gs は小さい軌道で通る。`qcod` の初期値・後退の違いの可能性があるが未解決。
- `Polarization` は k64 モードでは NaN (§6)。v2 モードの有限値の限界は §8a。`Polarization2/4/6` は実装していない。
- ROTATE 付き空洞の `Emittance[]` のスピン回転は照合していない (§5.3)。
- 追跡は (k64 モードでは) 損失の処理なし。numba の高速追跡 (`TrackRun`) を使わないので遅い (Python の追跡)。
- `tbrote` の更新 (ずれ・DROTATE・CHI2 のある偏向磁石) は gs で確認していない (B56)。
- `Emittance[]` が SOL を含むときの偏極の量は上記の SOL の例 (四極 1 つ + ドリフト) だけ照合している。SOL の中の BEND・MULT・CAVI がないリングに限る。
- スロー経路 (`PYSAD_NO_ELKERN=1` 相当) の `_tquade` は、ずれだけで回転のない四極のヤコビアンを、ずらした後の座標で評価する (ずらす前だと POL は numba の速い経路を使えないので εy に 7e-5 の誤差が出る)。

## 8a. v2 モードで変わった点 (k64 モードは上の §1〜§8 のとおり)

試験 `test_v2_A_spin.py`・`test_v2_A_tchge.py` (放射・エミッタンスの記録 A)、`test_v2_C2_spin.py` (スピン追跡の記録 C2)。記録は `docs/history/v2_changes_A.md`・`v2_changes_C2.md`。根拠は解析式・独立な数値計算 (Thomas-BMT 方程式の積分、モンテカルロ、FLUC つきのスピン追跡) で、gs ではない。

### Polarization (平衡偏極度) の有限値 (B55)

- **K64 の動作**: `spdepol` の `exp(cb²−cφ)·cerfc(cb)` が 0·∞ で NaN。桁あふれだけ直して有限にしても、誤差のあるリングで追跡と 3〜5 桁合わない (側帯の項に 2π が余分に掛かり、ほかにも 1/減衰率 の大きさの因子がある)。式の導出に戻らないと直せないので、v2 はこの式を使わない (`spin_results` の `equpol_k64` に桁あふれだけ直した値を参考として残した)。
- **v2 の動作**: (1) `exp(cb²)·erfc(cb)` は Faddeeva 関数 w(i·cb) で計算する。(2) 平衡偏極度は、線形のスピン軌道結合から脱偏極率 `rate` を求め、P = (8/5√3)·sdamp/(sdamp+rate) とする。`rate` は、要素ごとの軌道の行列 T_i、放射の拡散行列 B_i、スピンの回転 R_i、スピン軌道結合 D_i (Emittance[] の srot) から、スピンの周期解 n(u) = n0 + N u を N_{i+1} = (R_i N_i + [n0]× D_i) T_i⁻¹ で解き、rate = ½ Σ_i tr(N⊥ B_i N⊥ᵀ) (1 周あたり) とする。Derbenev–Kondratenko の式の (11/18)|γ∂n/∂γ|² の項と同じ (`pysad/spin.py::slim_depolarization`)。
- **結果**: 平面リングで **92.376% (Sokolov–Ternov)**。誤差のあるリング (12 セルのリングに四極の縦ずれ ±3 mm、閉軌道 y 6 mm) は、FLUC つきのスピン追跡 (RAD FLUC POL、16 粒子) で測った脱偏極率と下表のとおり合う:

| エネルギー・RF | スピンチューン (端数) | 追跡で測った脱偏極率 [/周] | v2 の `depol_rate` [/周] | v2/追跡 | K64 の式 (桁あふれだけ直す) |
|---|---|---|---|---|---|
| 3 GeV、1.5 MV | 0.403 | 3.4e-11 | 5.5e-11 | 1.6 | 2.2e-5 |
| 2.903 GeV、1.5 MV (νx = 0.290 の近く) | 0.293 | 1.05e-9 | 1.28e-9 | 1.2 | 1.7e-5 |
| 4.2 GeV、12 MV | 0.236 | 4.4e-8 | 2.3e-8 | 0.53 | 4.3e-5 |

  追跡の統計は ±25% (16 粒子) で、測った 1 − ⟨S·n⟩ は (4.1 ± 0.9)e-7、(8.4 ± 2.0)e-6、(1.31 ± 0.27)e-4。平衡分布での n の広がり `spin_spread` は 3 GeV で 6.4e-6 で、追跡の 1 − ⟨S·n0⟩ (FLUC 中、平衡後) 3e-6〜1.2e-5 と合う。線形モデルのモンテカルロでは、1 − ⟨S·n(u)⟩ が 20000 周まで式の rate (8.3e-11 /周) どおり直線で増える (±2%)。5 GeV・8 MV の例は追跡で粒子がバケットから落ちて比較にならなかった (4.2 GeV・12 MV に替えた)。
- 出力: `emit6d.emittance6d(lat, p, rfsw=True, radcod=True, pol=True)['spin']` に `equpol` (v2 の偏極度)、`depol_rate` (1 周あたりの脱偏極率)、`taudep` (脱偏極時間 [s])、`spin_spread`、`spin_orbit` (入口の ∂n/∂u、3x6)、`equpol_k64` (K64 の式の参考値) が入る。SADScript の `Polarization/.em` が v2 の値になる。
- **限界 (使えない範囲・未確認)**: 軌道の振幅について 1 次 (SLIM 相当)。追跡との一致は 3 例で 0.5〜1.6 倍 (16 粒子、統計 ±25%、確かめたのは小さな 1 つのリング)。シンクロトロン側帯などの高次の共鳴、運動論的偏極 (b·γ∂n/∂γ の項) は入れていない。エネルギー幅の大きいリング (aγσδ/νs が 1 程度以上) では側帯の寄与が増えるはず。K64 の spdepol (2 次の摂動、シンクロトロン側帯) は大きさが合わないので置き換えた (K64 の式を直すには導出が要る)。SOL 領域内の MULT・BEND の POL (Emittance[]) は未対応のまま。POL のときは四極・MULT を k·(L/n)² ≤ 2e-3 に分割するので Emittance[] が遅くなる。`Polarization2/4/6` は実装していない。

### Emittance[]・追跡のスピンの修正 (K64 の誤り)

| ID | K64 の動作 | v2 の動作 | 根拠 |
|---|---|---|---|
| B58 | BEND の端の縦磁場の積分 bsir0: 出口の項 (bsi2) が常に捨てられる、K1 の項に積分値 akx を掛ける、スピン軌道結合に端の項の y・x 依存が入らない。追跡の tbendi は y/ρ がなく k·x·y を全スライスで足す | 入口で +(y/ρ + k x y)、出口で −(y/ρ + k x y) (k = K1/L)。軌道についての微分もスピン軌道結合に入れる。追跡の tbendi も同じ | 1 要素のスピン回転が Emittance[] と追跡で一致: K1=0.3 の BEND 1.7e-6 (K64 2.0e-3)、K1=0 の BEND 7.5e-5 (K64 1.0e-3)。スピン軌道結合は追跡の差分と 2% (K64 は 1/6) |
| B64 | tmulke の bsir0 が imag(cx)/aln (次元 1/m)。四極 (tquade) は入れない。四極・MULT の分割は放射の大きさだけで決まる | ψ = Im Σ k_n z^(n+1)/(n+1)! を入口で足し出口で引く。四極は k x y。POL のときは四極・MULT を k·(L/n)² ≤ 2e-3 に分割 | MULT (K2=30) 3.1e-7 (K64 9.4e-5)、MULT (K1・K2・SK2) 3.8e-6 (3.1e-4)。四極の分割が 1 個だと歳差の微分が追跡と 50% 違い、細かくすると 2% |
| B56 | tbrote の srot の更新が行の逐次更新で行列積にならない | (I+rr)·srot | 直交性が 1e-14 で保たれ、入口と出口が互いに逆 (K64 は 1e-5 以上崩れる) |
| B39 | tradke (SOL 内) の dbeam(8) が B12 − B44 + b·B13 | B12 − B34 + b·B13 | 力学的運動量の共分散 |
| B81 | sprot の縦磁場 bsi: δ≠0 で (1+a)B∥ の歳差が 1/(1+δ) 倍にならない (δ = −5e-3、0.12 rad で 6e-4 rad)、軌道が傾くと bsi·ẑ の運動量に垂直な成分が B⊥ に二重に入る、双極 + ソレノイドは T-BMT と 3.6e-3 rad 合わない | SOL 領域と UND では、運動量の向きの磁場の積分 ∫B∥dl/(Bρ) を運動量の向きに置く | Thomas-BMT の積分 (v2 1e-7〜1e-5、K64 の式 6e-4〜4e-3) |
| B82 | 領域内の四極 (RAD): bsi を最初のスライスは半分、最後のスライスは半分で上書きするので 1 スライス分足りない。放射の基準の向きが正準運動量のまま | スライスごとに bz·aln/(1+δ)。放射の基準は力学的運動量 | 軸上の粒子のスピンの回転 (1+a) bz L: K64 は半分、v2 は 5e-15。軸外の放射損失: v2 の四極と MULT (K1) が 3e-4、K64 は 2.2 倍 |
| B83 | sprot の歳差の角が (1+aγ)·2 sin(θ/2) (刻みが粗いと θ²/24 小さい。0.26 rad の BEND を 3 刻みで 3e-4) | 厳密な θ。**すべての要素の追跡のスピンに効く** | 平面リングの 1 周の回転角 = 2π aγ (K64 の式は 4e-3 rad ずれる)。リング + ソレノイド (0.5 GeV、aγ = 1.135、1 m の 0・1・2.5 T): cos(μ/2) = cos(π aγ) cos(φ/2)、φ = (1+a) bz L と 2e-5 |

(B81〜B83 は v2 の作業で採番した K64 のバグ。docs/21)

### TrackParticles のスピン追跡 (track6 の経路。POL + RAD)

- 領域の外の要素: 上の B83 (刻みの角)、B58、B56 が入る。K1=0 の BEND の追跡のスピン (tbend + tradk) は Emittance[]・tbendi (v2) と 2.2e-5 rad 違う (放射の分割が粗く、スライスの前のエネルギーで歳差を計算する。K64 の精度の範囲で、直していない)。
- **損失**: POL + RAD の TrackParticles も tapert で落ちる (最初の要素・15 要素ごと・区間の終わり。LOSSMAP の 2 行も付く。発散は例外にせず落とす)。gs: 落ちた粒子の要素 74・周 1・座標 .022037636775、生存粒子のスピンが一致 (k64 モードで 1e-12。この項目は k64 モードにもある)。MAP 要素や空間電荷が落とした粒子の要素番号は 0。
- **SOL 領域** (v2 モード): 領域内の DRIFT・QUAD・BEND・MULT・CAVI・SOL (ずれ・回転・F1 つき) で通る。Thomas-BMT の積分 dS/dl = S × [(1+aγ)B⊥ + (1+a)B∥]/|p| と 1 GeV・0.4 m (軸上と軸外、δ = ±5e-3、3 つの向き) で比べた: DRIFT 5e-7、BEND 3e-6、BEND + bz 1e-5、QUAD + bz 3e-6、MULT (K0+K1) + bz 5e-5。スピンと運動量の角 = aγ × 偏向角 (2e-4)。座標は高速追跡と 1e-10。SOL の端の回転はスピンも回す (s ← R s)。歳差は放射の刻みごとに 1 回の回転にまとめるので、一様磁場の刻みは厳密、四極・多極の中は刻みの粗さの誤差 (`EPS` を小さくすると細かくなる)。k64 モードでは領域内は DRIFT・QUAD だけ。**領域内の APERT のスピン追跡は NotImplementedError。領域内の MULT・CAVI の POL は未確認。**
- **UND** (磁場つき): 刻みごとに [写像 → 放射とスピンの歳差]。縦磁場 Bs (軸外) と、周期の途中で終わる端のポテンシャルの飛び (薄い蹴り) の歳差も入る。B = rot A の中の Thomas-BMT の積分と 2e-6 (POLE = 0.25、5、5.25、中点面の外 2〜3 mm)。k64 モードは NotImplementedError。
- **PHSROT**: 位相空間の写像 (減衰・拡散) だけで、スピンは動かさない (磁場を持たないので歳差が決まらない)。リング全体を PHSROT 1 個で表すとスピンチューンは 0 になるので、偏極の計算には使えない。
- **Emittance[]・EMIT の POL の量 (SpinTune・Polarization など) は、SOL 領域の中は K64 の式 (B81 の bsi) のまま**。軸上のソレノイドでは追跡と 2e-5 で合うが、ずれて傾いたソレノイド (軌道が軸を外れる) ではスピンチューンで 2.3e-4、n0 の向きで 4e-3 rad 違う。直すには `spin.tradke_spin` の `xpzb = xpz + (bsir0 + bzh·2·al)·pr²` と 6 次元の微分 dxpzb を一緒に直す必要があり、未着手。領域の外の Emittance[] のスピンチューンは B83 の分だけ追跡 (v2) と違う (平面リングで 2.8e-5)。

## 9. 関連文書

- `docs/36_spin_csr.md`: POL と CSR の要約 (§1 とスピン追跡の節)
- `docs/29_sol_optics.md`: SOL の Twiss・Emittance・追跡
- `docs/07_emittance.md`、`docs/08_tracking.md`: Emittance と TrackParticles の一般説明 (POL の記述は更新が必要な箇所がある)
- `docs/21_sad_bug_report.md`: K64 の癖 B55〜B58
- `docs/history/work_unimpl.md`: POL・SOL・追跡の作業記録と未解決の差
- `docs/45_radlight_guide.md`: RADLIGHT (`TrackParticles` が軌道を返す)
