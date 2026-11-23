# v2.5.0 記録: tsolrot の傾いた枝の更新の符号 (B337) と MULT の CHI1・CHI2 の符号

日付: 2026-11-23。試験 `pysad/tests/script_tests/test_v250_tsolrot_tilt.py` (7 件)。gs の保存した出力 `pysad/tests/script_tests/scripts/v250_tsolrot_tilt__{plain,tilt}.*`、`v250_mult_chi_sign__{none,chi1p,chi1m,chi2p,chi2m}.*`。

## 1. (E) 症状

`tsol.f` の `tsolrot` は SOL 領域内の MULT・QUAD の取り付け誤差 (DX・DY・DZ・CHI1・CHI2・CHI3) の座標変換 (MULT の `tmulti` は −chi1、−chi2 を渡す)。Bz ≠ 0 かつ CHI1・CHI2・DZ ≠ 0 の枝 (`tsol.f:83-100`) は、傾いた一様磁場 |bz|·φ のらせん運動 u(θ) = u_l + u_t cos θ + u_b sin θ、r(θ) = r0 + [u_l θ + u_t sin θ + u_b (1 − cos θ)]/b_p を、面の条件で θ を求め (ニュートン法)、位置・運動量を更新する。`xsincos` の戻り値 dcosphi = cos θ − 1 を、面の条件 (`s = plz*xsinphi + pz0*sinphi - pbz*dcosphi`) では + pbz (1 − cos θ) の正しい符号で使うが、更新 (`x += (… + pbx*dcosphi)*alb`、`px -= ptx*dcosphi + pbx*sinphi …`) では 1 − cos θ の想定の符号で使う (古いコメントアウトの `dcosphi = 2 sin²(θ/2)` = 1 − cos θ の名残)。物理ガイドの確認 (numpy の写しとらせん運動の ODE) が指摘した不整合。

## 2. (E) 確認

- 基準 (独立): 長さ 0 の傾いた要素は恒等写像。入口と出口 (`tsolrot` の ent = false の枝) を numpy に写して合成する。更新を Fortran のまま (cos − 1): 恒等からのずれ 2.15e-7・3.87e-7・1.45e-6・6.48e-6 (3 GeV/c、bz = 1.3 /m、CHI1 = 0.03、CHI2 = −0.02、DZ = 0.01、4 粒子、px 最大 0.05)。更新の符号を反転 (1 − cos): 5e-18 以下。
- gs: SOL (BZ = 13 T、BOUND = 1) → MULT (L = 0、K2 = 1e-9) → SOL (BZ = 0)。傾きなしは入力に戻る (5e-16〜6e-15)。傾き (CHI1 = 0.03、CHI2 = −0.02、DZ = 0.01) は出力 − 入力が 2.19e-7・3.88e-7・1.43e-6・6.53e-6 で、numpy の「Fortran のまま」と 3% 以内で同じ大きさ。符号を反転した写しなら恒等になるので、gs は不整合のある式の挙動をしている (K64 の誤り)。CHI1 の符号を逆にしても同じ大きさ (2.15e-7〜6.48e-6)。
- pysad は SOL 領域内の MULT・QUAD の CHI1・CHI2・DZ の枝を持たない (未実装)。翻訳の誤りではなく K64 の誤り。

## 3. (F) MULT の CHI1・CHI2 の符号

gs の MULT (L = 1 m、K1 = 0.3、3 GeV/c) に CHI1・CHI2 = ±0.05 を付け、軸上の粒子 (x = y = px = py = 0) の出口を見る。

| 指定 | 出口 | gs | ヘルプの R₁・R₂ の軸の ODE |
|---|---|---|---|
| CHI1 = +0.05 | x、px | +1.19393e-3、−1.83961e-4 | +1.19541e-3、−1.83930e-4 |
| CHI2 = +0.05 | y、py | −1.30637e-3、−1.91410e-4 | −1.30801e-3、−1.91440e-4 |

CHI = −0.05 は符号が反転。ODE の基準: ヘルプの (8.10) の R = R3 R2 R1 で、要素座標 = R·(実験室座標 − d)。要素の軸 (要素座標の z 軸) の実験室での向きは R の 3 行目: CHI1 で (−sin χ1, 0, cos χ1)、CHI2 で (0, −sin χ2, cos χ2) (R2 = [[1,0,0],[0,c,s],[0,−s,c]])。中心 (s = L/2) を通るこの軸のまわりの厚い四極 (x 収束・y 発散、K = 0.3 /m²) の運動方程式を DOP853 で解いた。符号まで 0.13% 以内で一致 (残りは傾いた入口・出口の面の幾何の違い)。つまり K64 の MULT の CHI1・CHI2 の符号はヘルプの R₁・R₂ と同じ (`tmulti` の −chi は `tsolrot` の内部の規約)。参考: DX = +1e-3 は x = +1.46e-4 (要素が +x へずれる向きと整合)。

## 4. 限界・未確認

- (E) は 1 つの傾き・bz = 1.3 /m・領域の長さ 0。傾いた領域の中の L ≠ 0 の MULT の本体 (傾いた磁場の中の `tdrift`・`tsolqu` の bxs・bys) は見ていない。
- (F) は L = 1 m の四極 1 つ・CHI = ±0.05。DZ と CHI3 (ROTATE) の符号は見ていない。
- pysad の MULT の CHI1・CHI2・DZ は、領域の外 (Bz = 0) でも未実装 (指定しても結果が変わらない)。実装する場合: k64 は K64 のまま、v2 は B337 の符号を反転する。
