# v2 の変更 (担当 C): 空間電荷の追跡 (SPAC / PSPAC / SPCH / DRIFT の RADIUS)

2026-11-17。対象: `pysad/spch.py` (新規)、`pysad/track6n.py` と `pysad/ffscmd.py` のつなぎ、`pysad/parser.py` の 1 行、試験 `pysad/tests/test_v2_C_spch.py` (22 本)。
v1 までは「フラグを読むだけで、空間電荷は追跡に入らない」だった。v2 モードで、`TrackParticles` に空間電荷が入る。k64 モードは今までどおり入らない。

## 1. K64 の空間電荷のモデル (ソースを読んだ結果)

| フラグ | 入る場所 | K64 のソース | モデル |
|---|---|---|---|
| `SPAC` | `DRIFT` (`RADIUS>0`)。SOL 領域の中のドリフトも | tturn.f:288 → spkick.f `spdrift_free`、tsol.f:639 → `spdrift_solenoid` | ドリフトを 5 cm (`alstep`) ごとに割り、その間に、バンチの静止系 (ζ = γ z) での 3 次元の電位のキックを入れる。メッシュは粒子に合わせて張り、円筒の壁 (半径 RADIUS) は周上 12 点の像の電荷で近似する。`RADIUS<0` は半径 \|RADIUS\| の絞りだけ、`RADIUS=0` は何もしない |
| `PSPAC` | `SPCH` 要素 | tturn.f:521 → tspch_.c `tpspac` → spch.c / psn.c | 横 2 次元のポアソン方程式 (PSPACNX×PSPACNY の矩形メッシュ、境界で電位 0) の電位 Φ(x,y) と、縦の線密度 λ(z) の積。前の SPCH からの距離の分のキックを SPCH の位置で入れる |
| `WSPAC` | すべての要素の中点 | tintrb.f `twspac` | ビーム行列 (Emittance の結果) のガウス分布の力。**この項目の対象外 (pysad は未翻訳のまま)** |

`tspac.f` の `tspac` (粒子の rms からのガウス分布のキック) は、K64 のどこからも呼ばれない死んだコード (同じファイルの `bbkick1` だけが `GaussianCoulomb` と `twspac` 系で使われる)。

## 2. K64 の誤り (確かめた順)

| ID | K64 の動作 | 根拠 | v2 |
|---|---|---|---|
| SPC1 | **SPAC の力が 1/γ 倍に弱い**。spkick.f `spfield` は静止系の電位 φ' (= Σ 個数/r') の勾配をそのまま実験室系の場として使う (`ex = −∂φ'/∂x`、`ez = −∂φ'/∂ζ / γ`)。正しくは E⊥ = γ E'⊥、E_z = E'_ζ で、どちらも γ が 1 つ足りない | ソース + gs の実測 (下の表)。5 MeV で 1/9.8、20 MeV で 1/39 | 正しい係数 (§3) |
| SPC2 | SPAC のキックの後、`px(i)=px(i)/pr` と、すでに p で割ってある px をもう一度 (1+δ) で割る (spkick.f:199-200)。δ≠0 の粒子の横の運動量が 5 cm ごとに 1/(1+δ) 倍になる | ソース (spdeflect の入口の px は tturn の内部変数 px/p) | 正準運動量に戻して足し、新しい δ で割る |
| SPC3 | `MULT` の `RADIUS` の空間電荷は、`spkick` の呼び出しの引数が合っていない (tmulti.f:222、541 が `sx,sy,sz` を余分に渡す) | ソース | 未対応 (メッセージを出す)。§6 |
| SPC4 | SPAC のメッシュ: `nint(5.d0*dble(np/1000)**0.67d0)` が整数の割り算で、np<1000 では縦 10 枚に固定。横は `max(4, nint(npz^(1/4)/2))` で、np=2000 でも 8×8 | ソース。gs の結果は x と y が 16% 違う (下の表) | 24〜64 のノード、2 次の B スプライン |
| SPC5 | **PSPAC の既定 (`PSPACNZ=1`) では gs が落ちる** (rhocal で SIGSEGV)。縦のメッシュが 1 枚で、λ(z) の配列の外を読む | gs で確かめた (§5 のリングの例。`PSPACNZ=64` を与えると走る) | `PSPACNZ<4` はバンチに合わせた 64 ノードを自動で張る |
| SPC6 | PSPAC の λ(z) の配分 (spch.c `rhocal`): セルの前半では `rhoz[ii]+=rzi; rhoz[ii]+=rzin;` と両方を同じノードに足す (後半は `ii` と `ii+1`)。重みが z=セルの中央で不連続 | ソース | 2 次の B スプライン |
| SPC7 | PSPAC の周の終わりのキックは `pos(lend)-sspac0` で、`sspac0` は最後の要素の中点 (WSPAC 用の変数)。最後の SPCH から終点までの分が入らない。始めは `sspac1=0` で、区間を途中から始めても s=0 からの長さになる | ソース (tturn.f:201、270-277、584-597) | 区間の始めから数え、終わりに最後の SPCH からの残りを入れる |
| SPC8 | PSPAC の状態 (周の数、メッシュ) が C の static 変数で、`TrackParticles` を呼び直しても続きから数える。`PSPACNTURN` 周で解放した後も使う | ソース (tspch_.c) | 呼び出しごとに作り直す |
| SPC9 | `tspac` (死んだコード) の係数は `1/(p0 (1+p²))` = 1/(βγ·γ²) で、β が 1 つ足りない | ソース | 使わない |

PSPAC の式そのもの (spch.c `SpaceChargeBPSet` の K = 2 N r_0/(β²γ³)、`epotkick`) は正しい。v2 の PSPAC は同じモデルで、gs の結果とも合う (§5)。

**gs での SPAC の確認** (2026-10-19 ビルドの gs。1 m のドリフト `RADIUS=0.02`、半径 1 mm の一様な円・長さ 0.1 m の一様なバンチ、2000 粒子、初期の角度 0)。解析式は包絡方程式 a'' = K/a、K = 2 r_e λ/(β²γ³) の数値積分 (rms = a/2):

| 条件 | gs (K64) rms x / y [mm]、傾き px/x | 解析式 (K) | 解析式 (K/γ) | pysad v2 |
|---|---|---|---|---|
| 5 MeV (γ=9.84)、1e11 個 | 0.6226 / 0.7254、0.406 | 1.6387、1.150 | 0.6453、0.432 | 1.623 / 1.621 |
| 20 MeV (γ=39.2)、1e12 個 | 0.5084 / 0.5036、0.0390 | 0.7197、0.575 | 0.5060、0.0236 | 0.716 / 0.715 |

gs は K/γ の式に近く (x と y の差、メッシュの粗さの範囲)、K の式とは 5 MeV で傾きが 2.8 倍、20 MeV で 15 倍違う。

**例 `example/05_tracking/flag_SPAC` の記述について** (既存の文書は触っていない。親に報告): 「gs でも SPAC の効果が出ない」という結果は、`NOSPAC;` `SPAC;` を `Do[...]` の中に書いたためで、フラグが切り替わっていない (式の中ではフラグの命令にならない)。最上位に書けば gs でも効果は出る (上の表)。同じ例の (B) (RADIUS=−0.4 mm、SPAC、NORFSW、2 周) の生存数は、v2 の pysad で 84/200 になり gs と合う (v1 は 200)。(A) を最上位のフラグに直すと、この例は NORFSW なので v2 はエラーにする (`RFSW` にすれば走る)。

## 3. v2 の動作 (`pysad/spch.py`)

### SPAC (`Solver3D`、`kick3d`)
- バンチの静止系の座標: ζ = γ0 β0 (z/β_i) (z/β_i = −c·到着時刻。β0 γ0 はバンチの平均)。
- マクロ粒子を 2 次の B スプライン (TSC) でメッシュに載せる。メッシュは粒子に合わせる (中心 = 重心、半幅 = min(最大のずれ, 4σ)。外の粒子は端に寄せる)。
- ζ 方向は周期 FFT。横は、モードごとの自由空間のグリーン関数 2 K0(k d) (k=0 は −2 ln d) の、セルの上の平均 (原点の対数特異点は ∫K0(κu)u du = (1−κK1(κ))/κ² で厳密に積分) の畳み込み。
- 円筒の壁 (半径 R、軸は設計軌道、電位 0) は、像の電位 −2 Σ_m ε_m K_m(kR)/I_m(kR) · I_m(kr) I_m(kr') cos m(θ−θ') (m ≤ 40) を足す。k=0 は −(2/m)(r r'/R²)^m。
- 周期の像までの距離は min(3R, max(10×横の半幅, 3×縦の半幅)) (壁の中では場が exp(−2.4 Δζ/R) で消える)。
- 形状関数のぼかし (TSC の配分と補間、セルの平均) は 2 階差分で 1 次まで戻す。
- 粒子が自分で作った電位の分 (メッシュに対する位置で決まる見かけの自己力) を引く。残りは、他の粒子との対のエネルギーの勾配で、ハミルトン系のまま。
- キック (長さ L、r_c = q² r_e m_e/m、w = PBUNCH/NP0):
  - Δ(px/p0) = −C ∂φ'/∂x、Δ(py/p0) = −C ∂φ'/∂y、C = r_c w L/(β0² γ0²)
  - Δγ = −r_c w L ∂φ'/∂ζ。到着時刻は変えない (SAD の z = −v t に v'/v を掛ける)。
  - 一様な長い円柱では x'' = K x/a²、K = 2 r_c λ/(β²γ³)。電場の反発 (∝γ) と磁場の引力の打ち消し (1−β² = γ⁻²)、時間 L/(βc)、慣性 γ で 1/(β²γ³) になる。
- ドリフトの割り方は K64 と同じ: ndiv = max(1, nint(|L|/SPACSTEP))、SOL の中は nint(|bz L|/1.5) との大きい方。半分のドリフト → キック → … → 半分のドリフト。キックの前と終わりに、半径 |RADIUS| の円の絞り (外の粒子と NaN の粒子を落とす)。
- マクロ粒子 1 個の電荷は PBUNCH/NP0 (NP0 は `TrackParticles` に渡した粒子数。落ちた粒子の分だけ電荷が減る。K64 と同じ)。

### PSPAC (`Solver2D`、`kick2d`)
- K64 と同じモデル。Φ(x,y): 5 点差分のポアソン方程式を DST で解く (ノード i は x = (i − NX/2)·DX、境界のノードで 0)。λ(z): `PSPACNZ≥4` なら K64 と同じ固定メッシュ (ノード k は z = (k − NZ/2)·DZ)、それより小さければバンチに合わせた 64 ノード。
- Δ(px/p0) = −(r_c w L/(β²γ³)) λ(z) ∂Φ/∂x / N、Δγ = −(r_c w L/γ²) (dλ/dz) Φ / N (Φ = N 個の −2 ln d の和)。
- L は前の SPCH (区間の始め) からの距離。区間の終わりに、最後の SPCH からの残りを入れる。メッシュの外の粒子は蹴らない (K64 と同じ)。
- `PSPACNTURNCALC = n` (>0): n 周目までの電位を SPCH ごとに保存し、その後は使い回す。

### 設定 (FFS の変数)
| 変数 | 既定 | 意味 |
|---|---|---|
| `PBUNCH` | 1e10 | 1 バンチの実粒子数 |
| `SPACSTEP` | 0.05 | SPAC のキックの間隔 [m] (v2 で追加。K64 は 0.05 固定) |
| `SPACNX`、`SPACNY`、`SPACNZ` | 自動 | SPAC のメッシュのノード数 (v2 で追加)。自動は縦 clip(2·round(n^(1/3)), 24, 64)、横 clip(2·round(√(n/縦)), 24, 64) |
| `PSPACNX`、`PSPACNY`、`PSPACDX`、`PSPACDY` | 128、128、1e-3、1e-3 | PSPAC の横のメッシュ (セル数、幅 [m]) |
| `PSPACNZ`、`PSPACDZ` | 1、1e-3 | PSPAC の縦のメッシュ。NZ<4 は自動 |
| `PSPACNTURNCALC` | 0 | 電位を計算し直す周の数 (0 は毎回) |

### つなぎ
- `track6n.py`: `turn(..., spch=None)`。要素の分岐の先頭で `spch.handles(kind, 位置)` なら `spch.element(...)` (DRIFT と SOL 領域のドリフト)。周の始めに `spch.begin`、終わりに `spch.end`。`TrackRun.spch` (既定 None)。
- `ffscmd.py` (`TrackParticles`): `SPAC` か `PSPAC` が ON で v2 モードのとき `spch.make(self, tr, npart)`。
- `parser.py`: `DRIFT` の `RADIUS` を属性にする (`Element["RADIUS","D1"]` が gs と同じ値を返す。今までは 0)。モードによらない。

## 4. 根拠と試験 (`pysad/tests/test_v2_C_spch.py`、22 本、25 秒)

| 試験 | 比べる相手 (解析式・独立な計算) | 条件 | 一致 |
|---|---|---|---|
| `test_uniform_cylinder_perveance_gradient` | E'_x = 2λx/a²、軸上の電位 λ(1+2 ln(R/a))、R によらない | a=1 mm、R=2 cm と 3 mm、32768 粒子 | 勾配 5e-3、電位 2e-3 |
| `test_gaussian_matches_bassetti_erskine` | λ(z)×Bassetti-Erskine、中心の勾配 2λ/(σx(σx+σy)) | σx:σy=10:3、65536 粒子 | 全体 5e-3 (rms 1.5%)、中心 2% |
| `test_short_bunch_sphere_3d` | 球対称のガウス分布の E_r (3 次元) | σ=1 mm、R=10 cm | 5e-3 (rms 1.5%) |
| `test_point_charge_in_pipe_potential` | 円筒の中の点電荷 (ベッセル級数 + 周期の像) | σ=R/40、z=3.5σ (自由空間より 4.5% 小さい) | 3e-3 |
| `test_pipe_image_force_offset_beam` | 像の力 2λ(x0,y0)/(R²−x0²−y0²) | ずれ 0.4R〜0.75R、3 方向 | 1e-4 |
| `test_longitudinal_field_g_factor` | E_z = −λ'(z)(1+2 ln(R/a)) (軸上)、−λ'·2 ln(R/a) (r=a) | ガウスの λ、±2σ | 1% |
| `test_frozen_field_kick_is_symplectic[3]`、`[2]` | Jᵀ S J = S (差分のヤコビ行列、J−I の最大 >0.5) | 6 次元、2 点 | 1e-7 |
| `test_total_momentum_nearly_conserved_in_free_space` | Σ Δp = 0 | 壁なし | 力の絶対値の和の 1e-3 |
| `test_spac_kv_envelope_in_drift` | 包絡方程式 a''=K/a+ε²/a³ の数値積分 (a(1 m)=3.996 mm、NOSPAC は 2.236 mm) | 5 MeV、1e11、K=5.98e-6、ε=2e-6 | 2·rms 5e-3、傾き 2% |
| `test_spac_zero_charge_and_gamma_scaling` | PBUNCH=0 でドリフトと同じ (1e-15)、1 キックの傾き K L/a²、5 MeV と 50 MeV の比 1e3 | | 1%、3% |
| `test_spac_solenoid_brillouin_flow` | a'' = −k_L² a + K/a。Brillouin 流 (BZ=0.0816 T) で半径一定、0.7 倍で 1.81 mm (空間電荷なしは 0.14 mm) | SOL 領域 1 m | 1% |
| `test_ring_tune_shift_laslett[SPAC]`、`[PSPAC]` | 自前の 2×2 行列のリングに線形レンズ 2 r_e λ/(β²γ³)/(σ(σx+σy)) を入れたチューン (Laslett)。Δν = −0.00305 (SPAC: ドリフト 6.4 m)、−0.00381 (PSPAC: 全周 8 m) | 4 セルの FODO、5 MeV、ν=0.754、ε=1e-6、σz=5 cm、2.5e8 個、4096 粒子 | 追跡/モデル = 0.970、0.990 (SPAC)、0.980、1.004 (PSPAC)。5% で判定 |
| `test_pspac_spch_kick_uniform_beam` | 1 キックの傾き K L/a²、SPCH の位置での分け方、固定メッシュと自動 | | 1% |
| `test_pspac_longitudinal_kick` | Δγ = −r_e N L/γ² λ'(z)(1+2 ln(1.0787 b/a)) (正方形の箱の等価半径) | | 2% |
| `test_spac_k64_mode_unchanged` | k64 モードは NOSPAC とビット単位で同じ | | 厳密 |
| `test_spac_negative_radius_is_aperture_only` | RADIUS<0 は絞りだけ、`Element["RADIUS",..]=` で切り替え | | 厳密 |
| `test_spac_zero_bunch_length_error` | z がすべて同じとき・`NORFSW` のときのエラーの文。`RADIUS<0` は `NORFSW` でも走る | | |
| `test_example_transport_line_script`、`test_example_ring_script`、`test_element_radius_readable` | §5 のスクリプトが走り、値が解析式・gs の PSPAC と合う | | 2%、2e-3 |

注意 (試験を作るときに分かったこと): 小振幅の試験粒子を何個も同じ場所に置くと、試験粒子同士の力 (マクロ粒子 1 個の電荷) でチューンがずれる。チューンの試験は試験粒子を 1 個ずつ追う。強度を 10 倍にすると、追跡/モデルは 0.94〜0.96 に下がる (ビームの不整合と非線形で、線形レンズのモデルの側の限界)。

## 5. 使い方の例 (どちらも試験で走らせている。gs でも走る)

フラグは FFS の最上位に書く (`Do[...]` の中では効かない)。粒子には z の広がりが要る。`RFSW` を ON にする (`NORFSW` は区間ごとに z を 0 にするので線密度が決まらない。v2 は `NORFSW` で空間電荷を入れようとするとエラーにする。空洞が無ければ `RFSW` ON でも軌道は同じ。`RADIUS<0` の絞りだけなら `NORFSW` でも走る)。

### 輸送路 (SPAC)
```
MOMENTUM = 5MEV;
DRIFT D1 = (L=1.0 RADIUS=0.02);          ! RADIUS>0 のドリフトに空間電荷が入る (真空ダクトの半径 2 cm)
MARK  ORIGIN = ();
LINE TL = (ORIGIN D1);
FFS;
RFSW; NORAD; NOFLUC;
USE TL; INS; CALC;
$FORM = "12.6";
np = 2000;
PBUNCH = 1e11;                           ! 1 バンチの実粒子数 (マクロ粒子 1 個は PBUNCH/np 個)
fr[x_] := x - Floor[x];
a = 1e-3;                                ! 半径 1 mm の一様な円、長さ 0.1 m の一様なバンチ
rr = Table[a*Sqrt[fr[i*0.6180339887+0.5]], {i, np}];
th = Table[2*Pi*fr[i*0.7548776662], {i, np}];
x0 = rr*Cos[th]; y0 = rr*Sin[th];
z0 = Table[0.1*(fr[i*0.5698402910]-0.5), {i, np}];
zero = Table[0, {i, np}];
beam = {x0, zero, y0, zero, z0, zero, Table[1, {i, np}]};
rms[v_] := Sqrt[Plus@@(v^2)/Length[v]];
NOSPAC;
r = TrackParticles[{1, beam}][[2]];
Print["RES NOSPAC rms x=", rms[r[[1]]], " y=", rms[r[[3]]], " alive=", Plus@@r[[7]]];
SPAC;
r = TrackParticles[{1, beam}][[2]];
Print["RES SPAC   rms x=", rms[r[[1]]], " y=", rms[r[[3]]], " alive=", Plus@@r[[7]]];
```
結果: NOSPAC 0.5003 / 0.4997 mm。SPAC は pysad v2 1.623 / 1.621 mm (解析式 1.639 mm。バンチの端の粒子を含む)、gs 0.6226 / 0.7254 mm (SPC1)。0.7 秒。

### 小さなリング (SPAC と PSPAC)
```
MOMENTUM = 5MEV;
DRIFT D1 = (L=0.8 RADIUS=0.03);
QUAD  QF = (L=0.2 K1=1.2) QD = (L=0.2 K1=-1.2);
SPCH  SP = ();
MARK  ORIGIN = ();
LINE CELL = (QF D1 SP QD D1 SP);
LINE RING = (ORIGIN 4*CELL);
FFS;
RFSW; NORAD; NOFLUC;
USE RING; CELL; CALC;
$FORM = "12.6";
np = 1000;
PBUNCH = 2.5e9;
SPACSTEP = 0.1;                  ! SPAC のキックの間隔 [m] (既定 0.05)
PSPACNZ = 64; PSPACDZ = 0.008; PSPACNX = 64; PSPACNY = 64; PSPACDX = 4e-4; PSPACDY = 4e-4;
fr[x_] := x - Floor[x];
gs[k_, c_] := Table[Sqrt[-2*Log[fr[i*k+0.5]+1e-12]]*Cos[2*Pi*fr[i*c]], {i, np}];   ! ガウス分布 (乱数を使わない)
eps = 1e-6; bx = Twiss["BX", 1]; ax = Twiss["AX", 1]; by = Twiss["BY", 1]; ay = Twiss["AY", 1];
u1 = gs[0.6180339887, 0.7548776662]; u2 = gs[0.5698402910, 0.8191725134];
u3 = gs[0.6710436067, 0.4502835858]; u4 = gs[0.8566748839, 0.7338918566];
x0 = Sqrt[eps*bx]*u1; px0 = Sqrt[eps/bx]*(u2 - ax*u1);
y0 = Sqrt[eps*by]*u3; py0 = Sqrt[eps/by]*(u4 - ay*u3);
z0 = 0.05*gs[0.7071067812, 0.3247179572];
zero = Table[0, {i, np}];
beam = {x0, px0, y0, py0, z0, zero, Table[1, {i, np}]};
rms[v_] := Sqrt[Plus@@(v^2)/Length[v]];
NOSPAC; NOPSPAC;
r = TrackParticles[{1, beam}, 1, 1, 3][[2]];
Print["RES NOSPAC 3 turns rms x=", rms[r[[1]]], " y=", rms[r[[3]]], " alive=", Plus@@r[[7]]];
SPAC;
r = TrackParticles[{1, beam}, 1, 1, 3][[2]];
Print["RES SPAC   3 turns rms x=", rms[r[[1]]], " y=", rms[r[[3]]], " alive=", Plus@@r[[7]]];
NOSPAC; PSPAC;
r = TrackParticles[{1, beam}, 1, 1, 3][[2]];
Print["RES PSPAC  3 turns rms x=", rms[r[[1]]], " y=", rms[r[[3]]], " alive=", Plus@@r[[7]]];
PSPACNTURNCALC = 1;              ! 1 周目の電位を凍結して使い回す
r = TrackParticles[{1, beam}, 1, 1, 3][[2]];
Print["RES PSPAC frozen rms x=", rms[r[[1]]], " y=", rms[r[[3]]], " alive=", Plus@@r[[7]]];
```
3 周後の rms x / y [mm]:

| | pysad v2 | gs (K64) |
|---|---|---|
| NOSPAC | 1.750 / 1.024 | 1.750 / 1.024 |
| SPAC | 1.776 / 1.041 | 1.750 / 1.034 (SPC1: 弱い) |
| PSPAC | 1.777 / 1.042 | 1.776 / 1.042 |
| PSPAC (凍結) | 1.777 / 1.042 | 1.779 / 1.041 |

v2 の SPAC (3 次元、ドリフトの中) と PSPAC (2 次元、SPCH にまとめる) が合い、gs の PSPAC とも 4 桁目まで合う。pysad は 3 秒。gs は `PSPACNZ` を与えないと落ちる (SPC5)。

もう 1 つの gs との比較 (PSPAC、直線、5 MeV・1e11・半径 1 mm の一様なビーム、`D(0.5) SPCH D(0.5) SPCH` の 4 要素、`PSPACNZ=32 PSPACDZ=0.005 PSPACDX=1e-4`): gs は rms x 1.461 mm、y 1.259 mm (x と y が 16% 違う。原因は調べていない)、v2 は 1.228 / 1.226 mm、薄いキックの解析式 (中央部) は 1.248 mm。

## 6. 使える範囲と条件、使えない範囲

使える (確かめた範囲):
- `TrackParticles` の中の空間電荷。SPAC は `DRIFT` (`RADIUS>0`) と SOL 領域の中のドリフト、PSPAC は `SPCH` 要素。
- 電子 5 MeV〜50 MeV、1e8〜1e11 個/バンチ、ビーム半径 0.1〜2 mm、真空ダクト 3 mm〜3 cm、バンチの長さ (静止系の σ) がダクトの半径の 1/100〜50 倍。式は γ によらず正しい形 (静止系の静電場) なので、他のエネルギーでも係数は同じ。他の粒子 (`CHARGE`、`MASS`) は r_c = q² r_e m_e/m で入るが、電子以外は試験していない。
- 分布は任意 (一様、ガウス、楕円、ずれたビーム)。マクロ粒子は 1000 個以上を勧める (試験は 1000〜65536 個)。
- 精度の目安: 場は 0.5% (65536 粒子)、4096 粒子の既定のメッシュで中心の勾配 (チューンシフト) が 1〜3% 小さめ。
- 速さ: 1 キックが 9 ms (1000 粒子)、16 ms (4096)、180 ms (3.3 万)、0.7 秒 (13 万)。SPAC はドリフト 1 m あたり 20 キック (`SPACSTEP` で変えられる)。

条件:
- 粒子に z の広がりがあること (すべて同じ z ならエラー)。`RFSW` が ON であること (`NORFSW` はエラー)。z の広がりが極端に小さい (nm 程度) と、その線密度のとおりの巨大な力になる。
- SPAC の壁は、設計軌道を軸とする円筒。ビームの外縁が 0.9 R を超えると像の精度が落ちる (メッセージを出す)。
- PSPAC は「横の分布が z によらない長いバンチ」のモデル (バンチ長 ≫ ダクト/γ)。メッシュ (`PSPACNX·PSPACDX`) をビームに合わせて決めること (既定は ±64 mm、幅 1 mm。細いビームでは粗すぎる)。メッシュの外の粒子は蹴られない。壁は矩形 (メッシュの縁)。
- マクロ粒子の数が電荷を決める (PBUNCH/粒子数)。`TrackParticles` に渡す粒子の全部が 1 つのバンチ。

使えない (未対応。調べた内容):
- `MULT` の `RADIUS` (SPC3)、四極・偏向の中の空間電荷。SPAC はドリフトの中だけなので、リング全体の空間電荷は、ドリフトの占める割合しか入らない (試験のリングで 80%)。全周の分は PSPAC + SPCH で入れる。
- `WSPAC` (ビーム行列からの力、`SELFCOD`)、光学計算・`Emittance` への空間電荷。
- MAIN の `TRACK`、`DynamicApertureSurvey`、スピン (`POL`) の追跡、`TrackRun` で扱えない要素の並び (エラーにする)。
- k64 モード (空間電荷なし。§7)。
- 複数バンチ (`NBUNCH`)、壁の抵抗、円でないダクト (SPAC)。
- 未確認: `TRPT` の加速 (空洞の後の p0 を使うようにしてあるが試験していない)、ウェーク・放射との組み合わせ、SOL 領域の中の四極 (空間電荷はドリフトだけ)、バンチの静止系での速度が相対論的になる場合 (γ·x' ≳ 1 や大きな運動量の広がり。静電場の近似が崩れる)。
- gs との比較は §2 と §5 の 4 例だけ。「小さい RADIUS>0 で gs が落ちる」(docs/46) は再現を試していない。

## 7. k64 モードの扱い

今までどおり、SPAC・PSPAC は追跡に入らない (`test_spac_k64_mode_unchanged`)。理由: K64 の SPAC は SPC1・SPC2・SPC4 の誤りを持ち、それをビット単位で再現する翻訳 (粗いメッシュ、12 点の像、SOR の反復) は行っていない。PSPAC は式が正しく、v2 の値が gs と 3〜4 桁合うが、ビット単位ではない (配分の重み SPC6、双 3 次補間、FACR の解法が違う)。
`Element["RADIUS", ..]` が定義の値を返すようにした変更 (`parser.py`) はモードによらない (gs と同じになる)。

## 8. 要判断

1. k64 モードで PSPAC だけ v2 の実装を使うか (gs に 3〜4 桁で近づくが、「k64 = 空間電荷なし」という今の動作が変わる)。今は入れていない。
2. `example/05_tracking/flag_SPAC` と、docs/46・docs/01 などの「SPAC は読むだけ」「gs でも効果が確認できない」の記述の更新 (例は `Do` の中のフラグが原因。§2)。
3. `MULT` の `RADIUS` の空間電荷を v2 で入れるか (K64 は呼び出しが壊れている。入れるなら MULT の分割の間にキックを挟む)。
4. `WSPAC` (twspac) は式が正しい (係数 1/(β²γ³)) ので、ビーム行列を使う翻訳は別項目として可能。
