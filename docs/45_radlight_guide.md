# 45. RADLIGHT と放射場 (RadiationField・RadiationSpectrum) の解説

最終更新: 2026-11-17

粒子の軌道を記録して、観測点での放射場とスペクトルを計算する一連の機能 (RADLIGHT、`RadiationField`、`RadiationSpectrum`) と、関連する PHOTONS (光子のリスト) の解説。初めて読む人向けに、使い方、K64 の癖、gs との一致、限界をまとめる。
K64 のソースは `tlinit.f` (`tlinit`、`tlstore`、`tlresult`、`tlcalc`、`tlspect`)。短い記録は `docs/36_spin_csr.md` の RADLIGHT の節、gs で RadiationField が落ちる件は `docs/21_sad_bug_report.md` B49。

## 1. 何ができるか

流れは 3 段階である。

1. `RADLIGHT;` フラグを ON にして `TrackParticles` を呼ぶと、粒子の軌道 `{t, x, y, z}` (GEO 座標) が要素ごとに記録されて戻る。
2. `RadiationField[軌道, 観測点]` が、観測点での電場・磁場・ポインティングベクトルの時間変化を返す。
3. `RadiationSpectrum[{τ, f}, {λ1, λ2, dλ}]` が、場 f(τ) の cos・sin 積分 (スペクトルの係数) を返す。

| 機能 | 状況 | 実装 |
| --- | --- | --- |
| RADLIGHT フラグ、TrackParticles の戻り `{結果, 軌道}` | 対応。k64 モードは gs と差 0 (K64 の癖: 横方向なし・BEND の中なし)。v2 モードは 3 次元の軌道で BEND の中も記録 (§3a) | `ffscmd._RadLight`、`FFS.track_particles` 付近 |
| `RadiationField` | 対応 (gs がこのビルドで落ちるので gs 照合なし) | `radfield.radiation_field` |
| `RadiationSpectrum` | 対応 (gs と全桁一致。NPARA の fork 並列なし、逐次) | `radfield.radiation_spectrum` |
| PHOTONS と `PhotonList` | **BEAMBEAM の光子だけ** 対応 | `beambeam.py`、`ffscmd` |
| BEND・四極など他の要素の光子 | 未対応 (BEND・曲げのある MULT を含むと SADError) | |

## 2. RADLIGHT フラグと TrackParticles の戻り

`RADLIGHT;` を ON にすると、`TrackParticles` は RAD と TRPT を ON にして追跡し (`rad = RAD or RADLIGHT`、`trpt = TRPT or RADLIGHT`)、戻りが変わる。

- 通常の戻り: `{ls または ld, {x, px, y, py, z, dp, flag}}` (flag は最後の行)。
- RADLIGHT のときの戻り: **`{{ls または ld, {x, px, y, py, z, dp, flag}}, 軌道}`**。ここで `軌道` は粒子ごとの表で、`軌道[[粒子, 1]]` が時刻の表、`[[粒子, 2]]`・`[[粒子, 3]]`・`[[粒子, 4]]` が GEO 座標の x・y・z の表。
- 従って、1 粒子の軌道は `TrackParticles[d, 10][[2, 1]]` (SAD 付属の Oide の例の使い方) で `{{t...}, {x...}, {y...}, {z...}}` になり、そのまま `RadiationField` の第 1 引数にできる。結果だけがほしいときは `[[1]]`。
- 高速追跡 (`TrackRun`) で扱えない要素の並びでは `SADError` (`RADLIGHT: この要素の並びは高速追跡 (TrackRun) で扱えない`)。

### 動く SAD スクリプト例 (軌道の記録)

`pysad/tests/test_radlight.py` の台本 (gs と差 0。1 GeV の電子、3 セルのリング、2 粒子、3 周)。

```
MOMENTUM = 1e9;
DRIFT D1 = (L=1) D2 = (L=0.5);
QUAD QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
BEND B = (L=1.2 ANGLE=0.2617993877991494 E1=0.1308996938995747 E2=0.1308996938995747);
MARK ST = ();
LINE CELL1 = (QF D2 B D2 QD D1);
LINE RING = (ST 3*CELL1);
FFS;
USE RING;
CELL;
CALC;
RADLIGHT;
np=2;
x=Table[0,{7},{np}];
x[[1,2]]=1e-4;x[[3,2]]=-1e-4;x[[6,2]]=1e-3;
Do[x[[7,i]]=1,{i,np}];
d={3,x};
r=TrackParticles[d,9,1,2];
t=r[[2]];
Print["S1 ",t[[1,1]]];     (* 粒子 1 の時刻の表 *)
Print["S2 ",t[[1,2]]];     (* 粒子 1 の x の表 *)
Print["R ",r[[1,2]]];      (* 追跡の結果 *)
Exit;
```

### 動く SAD スクリプト例 (放射場とスペクトル)

SAD 付属の Oide の例の RADLIGHT の部分。偏向磁石を ANGLE=0 の「ステアリング」にして (K0=0.1)、BEND の中を記録しない K64 の仕様 (§3) を避けている。

```
MOMENTUM=1e9;
scale=1e-15;
tw=1e-16;
B ANGLE 0
B K0 0.1
FindPeak[l_]:=Position[Abs[l],Max[Abs[l]]][1,1];
RADLIGHT
np=1;
x=Table[0,{7},{np}];
Do[x[7,i]=1,{i,np}];
d={1,x};
r=TrackParticles[d,10][2,1];
f=RadiationField[r,{7.9918,0.249792,0}];
p=FindPeak[-f[3]];
f[1]=f[1]-f[1,p];
s=RadiationSpectrum[{f[1],f[3]},{0.001,1,0.001}*3e-7];
```

これは格子 (USE 済み) の続きで書かれたスニペットで、1 本で完結しない。この流れ (B を `ANGLE=0 K0=0.1` の STEER にして 1 周、RADLIGHT → `RadiationField` → `FindPeak` → `RadiationSpectrum`) は pysad で最後まで通る (2026-11-16。軌道 8 点 {4,8}、f が {10,8}、ピーク位置 3、スペクトル {3,1001})。gs は `RadiationField` で SIGSEGV (B49) するので通しの比較はできない。軌道 (TrackParticles の戻り) の部分だけ gs と一致する。複数区間の軌道は §3 の 5 の限界に当たる。

## 3. K64 の癖 (tlinit / tlstore / tlresult。k64 モードの動作。v2 モードは §3a)

1. **BEND の中を記録しない**。記録するのは「区間 (`tturn`) の始め」と「L≠0 で型が CAVI より小さい要素 (BEND 以外) の後」だけ (K64 `tturn.f:555`)。
   `BEND` と曲げのある `MULT` (ANGLE≠0) は後も記録されないので、偏向磁石の円弧に沿った放射が軌跡の点列に入らない。意味のある放射場を得るには、Oide の例のように B を `ANGLE 0` にするなどの工夫が要る。
2. **横方向の座標が軌道に入らない**。`tlinit` が `geo0` の軸 1・2 を設定しない (common ブロックの初期値 0 のまま)。軌道に現れるのは進行方向 (geo の第 3 軸) に沿った点だけ (gs で確認)。試験の台本では 2 行目 (軌道の x) に 0.3, 0.8, 1.3 … と進行方向の距離が並び、3・4 行目 (y, z) の表はすべて 0 だった。進行方向は局所の第 3 軸で、大域では x に出る (2026-11-16 確認: `geo[0]=[[0,0,1,0],[-1,0,0,0],[0,-1,0,0]]` で局所 z 軸の大域成分は (1,0,0)。進行方向の距離は軌道の 2 行目に出る)。geo の向きが別なら出る行も変わる。
3. 位置は straight (`dphi=0`) の積み上げで、BEND で向きが変わらない。
4. 時刻は `t = gt − z/(v0·(1−dv+dvfs))`。点は粒子ごとに保持され、損失した粒子の扱いは未確認 (損失のない試験だけ)。
5. 同じ位置の点が 2 つ並ぶことがある (区間の始めの点と、直前の要素の後の点。試験の S2 で 7.8 が 2 回出る)。**複数区間 (複数周) の軌道では区間の始めの点が重なり (t が同じ点が 2 つ)、pysad の `RadiationField` は `x must be strictly increasing sequence` で失敗する** (未修正の限界)。gs は `RadiationField` が SIGSEGV するので比較できない。
6. `RadiationField` は **gs (この K64 ビルド) で SIGSEGV する** (どの引数でも。`tlcalc` の `loc_sad(ktfaddr(ka)-3,dl)` が不正な番地を指すとみられる。B49)。そのため `RadiationField` の gs 参照は取れていない。

## 3a. v2 モードの軌道 (B65)

K64 の §3 の癖 1〜3 (BEND の中を記録しない、横方向の座標が入らない、BEND で向きが変わらない) を、v2 モードでは直した。試験 `test_v2_C_radlight.py` (3 件)、例 `example/v2_C/C5_radlight.sad`。

- **横の軸を幾何から取る**: 軌道の x・y・z は幾何 (`LINE["GX"/"GY"/"GZ"]` と同じ GEO 座標)。振幅 0 の粒子の点は幾何と一致する (ROTATE つき・負の角を含む。1e-11)。時刻 = 行路長/(βc) (1e-9)。
- **BEND は円弧で進め、中も記録する**: tlstore の dphi≠0 の式で進め、1/γ の角度あたり `RADLIGHTDIV` 点 (既定 16) を記録する。`RADLIGHTDIV = n;` で変える。高いエネルギーでは点数が γ に比例して増える (3 GeV・0.26 rad の BEND で 2.5 万点)。曲げのある MULT は出口だけ、CAVI も記録する。
- 中の点は、状態の写しを刻んだ BEND で進めて作る (本当の追跡は要素全体を 1 回)。DX・DY・F1・FRINGE のある BEND は中の点だけ 1e-6 m 程度ずれる。
- **根拠**: 偏向磁石の放射場のピーク E = e β²/(4πε0 ρ R (1−β)²) (Liénard-Wiechert) と `RadiationField` が 0.25% (RADLIGHTDIV = 16)、5e-5 (64)。条件は 10.2 MeV (γ = 20) の電子。**gs との比較は無い** (gs は `RadiationField` で落ち、k64 モードの軌道は別物)。
- **使えない・未確認**: SOL 領域を含む区間 (SADError。両モード)、COORD による幾何の変化、曲げのある MULT の中の点 (出口だけ)。複数周の軌道で区間の始めの点が重なり `RadiationField` が失敗する限界 (§3 の 5) が v2 でも残るかは再確認していない (要確認)。高エネルギー・多数粒子の規模は未確認。

## 4. RadiationField と RadiationSpectrum (`pysad/radfield.py`)

### RadiationField[traj, obs]

`traj = {{t1..tm},{x1..xm},{y1..ym},{z1..zm}}`、`obs = {x, y, z}` (観測点、GEO 座標)。Feynman-Heaviside の式の時間領域の場
`E = (μ0 e CHARGE / 4π)(c² n/R² + R/c · d(c² n/R²)/dt + d²n/dt²)`。時間微分は自然スプライン (K64 の `spline`、y''/6) で求める。
電荷は大域変数 `CHARGE` (既定 −1)。

戻り: `{{τ..},{Ex..},{Ey..},{Ez..},{Hx..},{Hy..},{Hz..},{Sx..},{Sy..},{Sz..}}` (10 本の表。τ は遅延を補正した観測時間、H は n×E/Z0、S は E×H)。
引数は 2 個、第 1 引数は 4 本のリスト、第 2 引数は 3 成分でなければ `General::narg` / `General::wrongleng` のエラー。

### RadiationSpectrum[{τ, f}, {λ1, λ2, dλ}]

場 f(τ) を区分線形近似して cos・sin 積分を解析的に足す (`tlspect`)。k = 2π/λ を λ2 から λ1 まで `na = int((λ2−λ1)/dλ + 2.5)` 点で振り、戻りは `{{k..},{c..},{s..}}`。
スペクトル強度は `c² + s²`。K64 の `NPARA` による fork 並列は未実装 (逐次で計算)。

## 5. gs との一致の度合い

| 項目 | 試験 | 条件 | 一致 |
| --- | --- | --- | --- |
| 軌道の記録 (RADLIGHT) | `pysad/tests/test_radlight.py` (参照は gs の出力) | 1 GeV、BEND を含む 3 セル、2 粒子、複数周、始点 ls≠1 | 22 点 × 4 本 (t, x, y, z) × 2 粒子で相対 1e-12 以下。実際は差 0 |
| `RadiationSpectrum` | `test_radfield.py::test_radiation_spectrum_matches_gs` (参照は gs の出力) | 60 点のガウス × cos、λ 1〜5 mm | 相対 2e-6 以内の試験 (gs と全桁一致) |
| `RadiationField` | `test_radfield.py::test_radiation_field_of_uniform_motion` | 等速直線運動 (β=0.8)、観測点 {0.3, 0.1, 2.0} | gs 参照なし。Liénard-Wiechert の速度場と内部の点で相対 6e-3 (スプライン微分の離散化誤差) |

注意: `RadiationField` は **gs と比べていない**。上の 6e-3 は物理の解析解との比較で、K64 の式 (微分の取り方) の翻訳が gs と同じ値を返すことは未確認。放射光 (BEND 内の円運動) の場、ピークの位置・強度なども未確認。

## 6. PHOTONS と光子

`PHOTONS;` を ON にすると、`TrackParticles` は光子のリストを作り、大域変数 `PhotonList` に入れる (K64 の `tphotonlist`)。要件と範囲は次のとおり。

- **RAD と FLUC が必要**。なければ `Track::photons:  PHTONS requires RAD & FULC in TrackParticles` (gs と同じ文面。綴りも同じ) で止まる。
- 光子の行は位置・方向・偏光・粒子番号・要素番号・エネルギー。
- **BEAMBEAM の光子 (ビームストラールング) だけ** 対応。`docs/24_beambeam.md` の試験で 7 個の光子の位置・方向・偏光・粒子番号・要素番号が一致、エネルギーは相対 3e-9。
- 他の要素 (四極・BEND・MULT など) が放射する光子 (`tphrec`) は含まれず、その乱数 (光子ごとに `tgauss` 2 回) も引かない。RAD FLUC で他の要素が光子を出すと、その後の乱数列が gs とずれる (3 周の台本で gs 43 個、pysad 12 個の例)。
  BEND と ANGLE≠0 の MULT を含む並びは SADError (`PHOTONS: BEND の光子 (tbend の放射) は未対応`)。
- RADLIGHT の軌道と PhotonList は独立の機能で、RADLIGHT の軌道は光子の位置とは無関係 (放射場 `RadiationField` は軌道そのものから場を作る)。

## 7. 限界・未確認

- **RadiationField の gs 参照がない** (gs は SIGSEGV で落ちる。B49。2026-11-16 にも引数を変えて再確認できず)。等速運動の解析解との比較だけ。K64 の式との照合は読みによるだけで、数値の一致は未確認。
- **k64 モードでは BEND の中を記録しない**ので、そのままの BEND では放射場の主要な部分 (偏向の円弧) が入らない (K64 の仕様どおり。v2 モードは §3a)。Oide の例の流れ (B を ANGLE=0、K0=0.1) は 1 周の軌道なら pysad で通る (§2)。複数区間の軌道では `RadiationField` が失敗する (§3 の 5)。gs は `RadiationField` で落ちるため通しの照合は無い。
- k64 モードでは横方向の座標が軌道に入らない (K64 の癖)。粒子の横位置のずれ (x、y) は軌道の x・y に反映されない (v2 モードは §3a)。
- 損失した粒子の記録、ずれのある要素、RF 加速 (TRPT) 中の軌道の記録は未確認。試験は損失なしの条件だけ。
- `RadiationSpectrum` の NPARA の fork 並列は未実装 (逐次)。
- 試験の条件はエネルギー 1 GeV (RADLIGHT)、観測点は 1 点のみ。他のエネルギー、近接場・遠方場の使い分けは未確認。
- PHOTONS は BEAMBEAM の光子だけ。BEND などの光子の出力と乱数列は未対応。

## 8. 関連文書

- `docs/36_spin_csr.md`: RADLIGHT の簡単な記録 (本書が詳細)。
- `docs/24_beambeam.md`: BEAMBEAM の光子 (PHOTONS)。
- `docs/08_tracking.md`: TrackParticles、`tradk` の放射。
- `docs/21_sad_bug_report.md`: B49 (RadiationField が gs で SIGSEGV)。
- `docs/history/work_unimpl.md`: RADLIGHT の作業記録。
