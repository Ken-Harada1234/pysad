# v2.0.0 の変更記録 (担当 A: K64 のバグと物理の誤り — 放射・エミッタンス・スピン・RADINT・tchge)

作成: 2026-11-17。対象: `pysad/radint.py`、`emit6d.py`、`elkern.py`、`perfkern.py`、`turnkern.py`、`lattice.py`、`calc6d.py`、`touschek.py`、
`emit_extra.py`、`spin.py`、`geo.py`、`parser.py`、`solmap.py`、`dapert.py`、`track6.py`・`track6n.py`・`track6n_ext.py` の放射・スピン部、`ffscmd.py` の該当箇所。

- 切り替えは `pysad/compat.py` (`PYSAD_COMPAT=v2|k64`)。**k64 では K64 (gs) と同じ動作** (既存の試験・回帰は k64 で走る)。
- v2 の値の根拠は、解析式・独立な数値計算・保存則・別のコードとの一致 (gs は根拠にしていない)。試験は `pysad/tests/test_v2_A_*.py` (67 件)。
- 使った条件: 3 GeV の FODO リング (24 セル 105.6 m / 91.2 m、12 セル 45.6 m / 62.4 m)、偏向磁石 L=1.2 m・15 度 (ρ=4.58 m)、RF 500 MHz。
  個別の条件は各試験ファイルの先頭に書いた。大きなリング・挿入光源・強い結合では確かめていない。

## 1. 一覧

| ID | K64 の動作 | v2 の動作 | 根拠 (式・数値) | 試験 | k64 モード |
|---|---|---|---|---|---|
| B06 | RADINT が全 BEND に最後の BEND の L/ANGLE を使う | その BEND 自身の L/ANGLE | 曲率の違う 2 種の BEND のリングで、自前の 3x3 行列の積分と I1〜I5 が 2e-5 で一致 (K64 は I1・I4・I5 が 1e-3 以上ずれる) | `test_v2_A_radint.py::test_per_bend_radius_v2` | 同じ (`stale_rho`) |
| B07 | RADINT の六極が直前の DL を使う。L=0 は 0.2 m。先頭が SEXT だと DL=0 で NaN。DKL2 を代入しない | 六極は自分の L と K2 を使う。L=0 は寄与なし。モジュール変数を毎回初期化 | 軌道 0 なら六極は I1〜I5 を変えない (K2=3 と 0 で 12 桁同じ)。NaN にならない | `…::test_sextupole_first_and_thin` | 同じ |
| B12 | RADINT の定数が単精度 (1e27、1.E-6、0.2、2.3026) | 倍精度 | ELOSS = Cγ E⁴ I2/(2π) と 1e-13 で一致 (K64 は相対 1.15e-8 小さい) | `…::test_double_precision_constants` | 同じ |
| B13 | RADINT の CAVI は最後の空洞の値 | VOLT は総和。FREQ・HARM は電圧のある最初の空洞 (違う周波数が混ざると注意を出す) | 2.5 MV + 1.5 MV → 4000 kV。νs = √(hαV cosφ/2πE) と 2e-3 | `…::test_cavity_voltage_is_summed` | 同じ |
| R1 (新) | RADINT の `QUANTUM LIFETIME(hours)` が exp(2.3026·τ[秒]) | τ[秒]/3600 | τq = τE/2·exp(ξ)/ξ、ξ = (δmax/σδ)²/2 | `…::test_quantum_lifetime_hours` | 同じ |
| R2 (新) | RADINT の `LINEAR CHROMATICITY` が −I9/4π (磁石の中だけの ∮γ ds: 色収差にならない。例: −2.36、実際は −7.30) | 運動量を ±1e-4 ずらした光学のチューンの差分 dν/dδ (六極・偏向の寄与を含む) | 同じスクリプトの `DP0` をずらした CALC の差分と 2e-5 で一致 | `…::test_chromaticity_is_tune_derivative` | 同じ |
| R3 (新) | RADINT は MULT・SOL を黙って無視する | 場のある MULT・SOL があれば「入っていない」と表示する (積分には入れない) | — | `…::test_mult_warning` | 表示なし |
| JX1 (新、docs/46 の「RADINT と EMIT の Jx の差」) | Emittance[] (tradke) と追跡 (tradk) は、偏向磁石の端の角 e の蹴りを偏向角 θ に入れるが、長さ al1 を変えない。放射 ∝ θ²/al1 なので I4 の端の項 −η tan(e)/ρ² が 2 倍になる。矩形磁石のリングで Jx = 1.097 (正しくは 1.006) | 端の角のある面では、場の中の行路 al1 を x·tan(e) だけ短くする (その微分も D 行列に入れる) | (a) 解析式 Jx = 1 − I4/I2 (自前の積分) と (b) 硬い端の厳密な幾何で運動量のずれた閉軌道を求め、場の中の行路 l(δ) から Je = 2 + dln l/dδ (Robinson) が 2e-6 で一致: E1=E2=0 で 0.92082、0.1309 で 0.94147、0.5 で 1.00575。**RADINT (0.941466) が正しく、K64 の EMIT (0.9651) が誤り**。v2 の Emittance[] は 1.00586 (EMITDIV=20 で 1.00580)、追跡の一周写像の固有値は 1.0033 | `test_v2_A_damping.py` (7 件) | 同じ |
| JX2 (新) | 追跡の tbendi (K1≠0 の BEND) は、放射で運動量が変わっても曲率半径 rhoe を入口の値のままにする。減衰の和 Jx+Jy+Je が 4 にならない (K1=0.1 で 4.05、Jx=0.7385) | スライスごとの放射のあとに rhoe = rhob·(1+δ) を更新 | 和が 4.000、Jx = 0.6875 (解析 0.68762、RADINT 0.68762、Emittance[] 0.6883) | `…::test_track_edge_and_sum_rule_v2` | 同じ |
| ALPHA (新) | Emittance[] の MomentumCompaction: RFSW OFF は −R56/C だけ (入口の分散の項 Σ R5j η_j が抜け、始点で値が変わる。FODO の例で 0.022898、正しくは 0.022137: 3.4% 大きい)。RFSW ON はシンクロトロンチューンから μ² で求める (μ²/12 だけ大きい: 0.25%) | RFSW OFF は −(R56 + Σ R5j η_j)/C。RFSW ON は 4 sin²(μ/2) | 自前の I1/C と 2e-4 (FODO の例は 6 桁: 0.0221368)。RFSW OFF のバケット高さ・バンチ長・νs (√α に比例) も直る: BucketHeight 0.009877 → 0.010045 (RADINT 0.010045、RFSW ON 0.010047)。docs/46 の「BUCKET HIGHT の差」はこの α の差 | `test_v2_A_radint.py::test_momentum_compaction_is_i1_over_c_v2` | 同じ |
| EMITDIV (pysad の未対応) | 大域変数 EMITDIV・EMITDIVB・EMITDIVQ (放射の分割の倍率) が Emittance[] に効く | 同じ (ゲートなし)。pysad は無視していた。MAIN 部の `EMITDIV = 20;` も効く | gs: `EMITDIV=20` で DampingRate が −.000286831320204 (pysad k64 と 12 桁一致) | `…::test_emitdiv_is_used`、`…::test_emitdiv_in_main_section` | 効く |
| B09 | tchge の入口で、ロール θ≠0 かつ ずれ≠0 (または DZ・DROTATE≠0) のとき、tbrote が回転行列を上書きし、一周行列に回転が入らない | 入口の回転を行列に掛ける (行列 = 軌道の写像の微分) | `_tchge_jac` が `_tchge_cod` の数値微分と 1e-8 で一致。ロール + ずれのある六極・BEND を含むリングの CALC のチューンが、追跡の一周写像の固有値と一致 (六極 1e-9、BEND 4e-6。K64 は 0.02〜0.18 ずれる) | `test_v2_A_tchge.py` (5 件) | 同じ |
| B56 | tbrote の srot (スピン) の更新が行の逐次更新で、行列積にならない | (I+rr)·srot | I+rr は回転行列 (直交、det=1)。直交性が 1e-14 で保たれ、入口と出口が互いに逆 (K64 は 1e-5 以上崩れる) | `test_v2_A_tchge.py` (3 件) | 同じ |
| B11 | 4D の CALC が縦方向の Twiss 列 (AZ BZ NZ ZX ZPX ZY ZPY) を作り直さず、6D の値が残る | 既定値 (AZ=0、BZ=1、他 0) | 4D では縦方向を計算しない | `test_v2_A_ffs.py::test_calc4d_resets_longitudinal_columns_v2` | 同じ |
| B15 | `UnstableLevel` が常に 0 | 不安定な点の数 nstab | MatchingResidual = 10·nstab + NetResidual (例: 2 と 20) | `…::test_unstable_level_v2` | 0 |
| B16 | CalculateOptics の stabz・tracez が未初期化 | idim=3 は tracez = 2cos(μz) と stabz、idim=2 は stabz=1・tracez=0 | Emittance[] のシンクロトロンチューンと 2e-3 で一致 | `…::test_calculate_optics_idim3_returns_z_trace_v2` | 1 と 0 (不定の代わり) |
| B14 | tsolqum (放射なしの MULT) の z の dv の項が −0.25·(ap−dv)·aln0 (= +0.25·dv·aln0) | tsolqu と同じ −dv·aln0 | 同じ K1・L の QUAD と MULT の z が 1e-12 で一致、K1→0 の MULT がドリフトと 1e-13 で一致 (K64 は 1.25·dv·L ずれる。20 MeV・L=0.5 m で 1.9e-5 m) | `test_v2_A_track.py::test_mult_z_equals_quad_and_drift_v2` | 同じ |
| B08 | tsolqur (RAD、SOL 領域内の四極) が最後の半ステップで最後の粒子の bzp を全粒子に使う | 粒子ごとの BZ/(1+δ) | 1 粒子の結果が、一緒に追う粒子に依らない (差 0。K64 は 2.6e-6)。放射なしとの差は 1e-8 以下 | `…::test_solenoid_quad_rad_is_independent_of_other_particles_v2` | 同じ |
| B39 | tradke (SOL 内) の dbeam(8) が B12 − B44 + b·B13 | B12 − B34 + b·B13 | 力学的運動量 (px + b y, py − b x) の共分散 A B Aᵀ の px・py 成分 | `test_v2_A_spin.py::test_kinetic_momentum_covariance_formula`、`…::test_tradke_px_py_diffusion_uses_b_y_py_v2` | 同じ |
| B58 | Emittance[] (POL) の BEND の端の縦磁場の積分 bsir0: 出口の項 (bsi2) が常に捨てられる。K1 の項は積分値 akx を掛ける (次元が合わない)。bsir0 を定数として扱い、スピン軌道結合 (srot(:,4:9)) に端の項の y・x 依存が入らない。追跡の tbendi は y/ρ がなく k·x·y を全スライスで足す | 入口で +(y/ρ + k x y)、出口で −(y/ρ + k x y) (k = K1/L)。その軌道についての微分もスピン軌道結合に入れる。追跡の tbendi も同じ端の項にした | ∫Bs ds/(Bρ) = ψ(x,y) (横磁場のスカラーポテンシャル)。軌道のずれた 1 要素のスピン回転が Emittance[] と追跡で一致: K1=0.3 の BEND 1.7e-6 (K64 2.0e-3)、K1=0 の BEND 7.5e-5 (K64 1.0e-3)。K1→0 で K1=0 に連続。スピン軌道結合 (z まわりの歳差の y 依存) は追跡の差分と 2% (K64 は 1/6) | `test_v2_A_spin.py::test_element_spin_rotation_emit_equals_tracking_v2`、`…::test_bend_spin_is_continuous_in_k1_v2`、`…::test_bend_spin_orbit_coupling_includes_fringe_v2` | 同じ |
| B64 | tmulke (Emittance[] の MULT) の bsir0 が imag(cx)/aln (横磁場そのもの、次元 1/m)。四極 (tquade) は入れない。四極・MULT の分割は放射の大きさだけで決まる (軌道が小さいと 1 個) | ψ = Im Σ k_n z^(n+1)/(n+1)! を入口で足し、出口で引く。四極は k x y。POL のときは四極・MULT を k·(L/n)² ≤ 2e-3 になる数に分割する (歳差の微分の精度のため) | 同上。MULT (K2=30) 3.1e-7 (K64 9.4e-5)、MULT (K1・K2・SK2) 3.8e-6 (3.1e-4)。四極の縦方向の歳差の微分は、分割 1 個だと追跡の差分と 50% 違い、細かくすると 2% で合う | 同上 | 同じ |
| B55 | `Polarization` が NaN (spintrm の exp(cb²−cφ)·cerfc(cb) が 0·∞) | (1) exp(cb²)·erfc(cb) を Faddeeva 関数 w(i·cb) で計算。(2) 平衡偏極度は、線形のスピン軌道結合 (DK の (11/18)\|γ∂n/∂γ\|² の項) の脱偏極率 rate から P = (8/5√3)·sdamp/(sdamp+rate) | 平面リングで 92.376% (Sokolov–Ternov)。誤差のあるリングは FLUC つきのスピン追跡で測った脱偏極率と 0.5〜1.6 倍で合う (§3。K64 の式は 3 桁以上大きい) | `test_v2_A_spin.py` (B55 の 5 件) | NaN |
| B63 | SOL (BOUND≠0) の長さが周長に入らない | 変更なし (**バグではない**)。pysad の S 座標のほうを直した | SOL に L のキーワードは無い (SADHelp)。L を書いても軌道・光学に入らない (追跡で x が動かない、チューンが同じ)。gs の LENG も L を入れない (46.6)。pysad は S 座標 (LENG) にだけ入れていた (46.9) | `test_v2_A_ffs.py::test_sol_has_no_length` | 同じ (ゲートなし) |
| B36 | L≠0 の TCAVI のリングの Emittance[] が gs で NaN | 変更なし (pysad は有限) | K0=0 で同じ長さのドリフトと 12 桁一致、K0 について連続 | `…::test_tcavi_with_length_emittance_is_finite_and_has_drift_limit` | 同じ |
| B21 | RFSW OFF・RADCOD OFF の閉軌道の δ | 変更なし | 「z が閉じる δ」(行路 = 設計周長) を 1e-15 まで解いている。追跡で 1 周して z のずれ 2e-10 m | `…::test_no_rf_closed_orbit_closes_z` | 同じ |

## 2. Jx の 2.4% の差 (docs/46 §4) の結論

- **物理として正しいのは RADINT (INTGRL.f) の Jx = 1 − I4/I2 (端の項 −η tan(e)/ρ²、係数 1)。** K64 の EMIT・Emittance[] と追跡は端の項が 2 倍になる。
- 理由: 端の角 e のある面では、外側 (x>0) の粒子は場に入るのが x·tan(e) だけ遅れ、場の中の行路がその分短い (磁場の強さは変わらない)。
  放射は B²·(行路) に比例するので、端の寄与は −x·tan(e)/ρ² (係数 1)。実の勾配 (四極成分) の 2kη/ρ の「2」は B² の変化から来るもので、端には付かない。
  矩形磁石 (e = θ/2 を両端) では、行路が x に依らず、I4 への偏向磁石の寄与がほぼ 0 になる。
- K64 の tradke は、区間の前後の運動量の向きの差 θ と長さ al1 = al − Δz から放射を θ²/al1 に比例させる。端の蹴り x·tan(e)/ρ は θ に入るが al1 は変わらず、
  θ²/al1 ≈ (al/ρ²)(1 − 2x·tan(e)/al) になる (正しくは 1 − x·tan(e)/al)。
- 独立な確認 (試験 `test_v2_A_damping.py` に入っている):

| 端の角の比 E1=E2 | 解析式 (自前の I4/I2) | 厳密な幾何 + Robinson | RADINT | Emittance[] k64 | Emittance[] v2 | v2 (EMITDIV=100) | 追跡 k64 | 追跡 v2 |
|---|---|---|---|---|---|---|---|---|
| 0 (扇形) | 0.920824 | 0.920824 | 0.920824 | 0.92234 | 0.92234 | 0.92092 | 0.91871 | 0.91871 |
| 0.1309 | 0.941466 | 0.941466 | 0.941466 | 0.9648 | 0.94263 | 0.94154 | — | — |
| 0.5 (矩形) | 1.005751 | 1.005751 | 1.005751 | 1.09735 | 1.00586 | 1.00580 | 1.0964 | 1.0033 |

- gs (K64、2026-10-19 ビルド) で同じリング (E1=E2=0.5) を走らせると、Emittance[] の DampingRate は {−.00028590723365, −.00026053029608, −.00049548682084}
  (Jx = 1.0974、pysad の k64 と 12 桁一致)、同じ gs の RADINT は DU = 1.00575100。K64 の中で EMIT と RADINT が食い違っていて、正しいのは RADINT。
  台本は `local/wake/A_jx.sad` (リポジトリ外)。
- 残る差: 既定の分割 (EMITDIV=1) では Emittance[] の Jx が 1.5e-3 大きい (放射の刻みが粗い: K64 の精度)。`EMITDIV=20;` で 1e-4 まで合う。
  Jx+Jy+Je は 3.9990 (4 との差 2.5e-4 = U0/2E の 2 次の効果)。
- 影響: K64 の Emittance[] は、端の角のある偏向磁石のリングで εx を (矩形磁石の例で) 8% 小さく、σδ を 2% 大きく出す。v2 は RADINT と合う。

## 3. 平衡偏極度 (B55) の確認

- K64 の spdepol の式は、桁あふれだけ直すと有限になる (`spin_results` の `equpol_k64` に残した) が、誤差のあるリングで追跡と 3〜5 桁合わない
  (下の表の「K64 の式」)。側帯の項に 2π が余分に掛かり (幅 → 0 の極限が幅 0 の式の 2π 倍になる)、ほかにも 1/減衰率 の大きさの因子がある。式の導出に戻らないと直せないので、v2 は使っていない。
- v2 の式: 要素ごとに、軌道の行列 T_i、放射の拡散行列 B_i、スピンの回転 R_i、スピン軌道結合 D_i (Emittance[] の srot) を取り、スピンの周期解 n(u) = n0 + N u を
  N_{i+1} = (R_i N_i + [n0]× D_i) T_i⁻¹ から解く。脱偏極率は rate = ½ Σ_i tr(N⊥ B_i N⊥ᵀ) (1 周あたり)。これは Derbenev–Kondratenko の式の (11/18)|γ∂n/∂γ|² の項と同じ
  (B_i の δδ 成分 = (55/24√3) r_e ƛ γ⁵/|ρ|³ ds、Sokolov–Ternov の率 = (5√3/8) r_e ƛ γ⁵/|ρ|³ で、比が 11/18)。`pysad/spin.py::slim_depolarization`。
- 確認 1: 要素ごとの写像 s' = (1 − [D u]×) R s が、追跡 (track6 のスピン) の 1 周のスピンの軌道依存を 2〜6% で再現する (`test_spin_orbit_response_matches_tracking_v2`)。
- 確認 2: 線形モデル (同じ T_i・B_i・R_i・D_i、スピンは厳密な回転) のモンテカルロで、1 − ⟨S·n(u)⟩ が 20000 周まで式の rate (8.3e-11 /周) どおり直線で増える (式の検算、±2%)。
- 確認 3: FLUC つきのスピン追跡 (RAD FLUC POL、16 粒子、閉軌道上・スピン n0 から開始) を n1 周行い、そのあと FLUC を切って軌道を閉軌道に減衰させ、残った 1 − ⟨S·n⟩ (不可逆な脱偏極) を測った。
  リングは `test_v2_A_spin.py` の 12 セルリングに四極の縦ずれ ±3 mm (閉軌道 y 6 mm)。台本と出力はリポジトリ外 (作業用)。

| エネルギー・RF | スピンチューン (端数) | n1 [周] | 追跡で測った 1−⟨S·n⟩ | 追跡の脱偏極率 [/周] | v2 の `depol_rate` | v2/追跡 | K64 の式 (桁あふれだけ直す) |
|---|---|---|---|---|---|---|---|
| 3 GeV、1.5 MV | 0.403 | 12000 | (4.1 ± 0.9)e-7 | 3.4e-11 | 5.5e-11 | 1.6 | 2.2e-5 |
| 2.903 GeV、1.5 MV (νx = 0.290 の近く) | 0.293 | 8000 | (8.4 ± 2.0)e-6 | 1.05e-9 | 1.28e-9 | 1.2 | 1.7e-5 |
| 4.2 GeV、12 MV | 0.236 | 3000 | (1.31 ± 0.27)e-4 | 4.4e-8 | 2.3e-8 | 0.53 | 4.3e-5 |

- 平衡分布での n の広がり `spin_spread` = ½⟨|N u|²⟩ は 3 GeV で 6.4e-6。追跡の 1 − ⟨S·n0⟩ (FLUC 中、平衡後) は 3e-6〜1.2e-5 で合う。
- 途中でわかったこと: (a) Emittance[] のスピン軌道結合に端の縦磁場の微分が無いと、BEND の y 依存が追跡の 1/6 になる (B58 で直した)。
  (b) 四極の分割が 1 個だと歳差の微分が粗く、脱偏極率が 1.5〜3.4 倍大きく出る (B64 の分割で直した。`EMITDIV` を上げても同じ効果)。
  (c) 5 GeV・8 MV の例は追跡で粒子がバケットから落ちて比較にならなかった (4.2 GeV・12 MV に替えた)。

## 4. 「要判断」(直していない。根拠が足りないか、仕様の判断が要る)

| 項目 | 見つけたこと | 直さなかった理由 |
|---|---|---|
| tradke の dbeam の b = bzh·0.5 (photons.f) | 力学的運動量は px + bzh·y なので b = bzh のはず (B11 + b(2B13 + b B33) の b が半分) | bzh の定義の別の約束がある可能性を排除できない。SOL 内の放射のビーム行列の「ビームの広がり」の項だけに効き、ε は全桁変わらない |
| tbgrote の trans(5,6) (tchge.f) | 写像の微分 zf/pzf·(1 − pr/pzf·r33·dpzidp) (コメントアウトされた式) ではなく r33·(dpzf·pr/pzf + dpzi/pzi) − dr33。dr=0 でも dpzi·(pr+1)/pzi が残る | 作者が意図して置き換えた式。2 次の小さい量 |
| RADINT の CKL1Y = −cos(θ)·DKL1 (BMAG)、−DKL1 (QUADRU) | sin の誤記に見える (結合のあるときの I4V に効く) | 結合のあるリングの独立な I4V の計算を持っていない |
| RADINT の BEND の K0 | 無視される (ANGLE だけ) | 仕様かもしれない |
| RADINT の BUCKET HIGHT と Emittance[] の BucketHeight (0.010045 と 0.010035) | 解決 (一覧の ALPHA)。バケット高さ ∝ 1/√α で、K64 の RFSW ON の α が μ²/12 = 0.25% 大きかった。v2 は 0.010047 と 0.010045 (残りの 2e-4 は空洞と放射の分) | — |
| K64 の spdepol (2 次の摂動、シンクロトロン側帯) | §3 のとおり大きさが合わない。v2 は線形の式に置き換えた | K64 の式を直すには導出が要る |
| v2 の Polarization の限界 | 軌道の振幅について 1 次 (SLIM 相当)。追跡との一致は 3 例で 0.5〜1.6 倍 (16 粒子、統計 ±25%)。シンクロトロン側帯などの高次の共鳴、運動論的偏極 (b·γ∂n/∂γ の項) は入れていない。エネルギー幅の大きいリング (aγσδ/νs が 1 程度以上) では側帯の寄与が増えるはず | 高次の共鳴は線形理論の外。確かめたのは小さな 1 つのリングだけ |
| CALC (tbende) と追跡 (tbend) の、ずれた BEND | 矩形磁石を DX=2 mm ずらすと CALC のチューンが 2.9e-5 変わる (追跡は変わらない。平行な面の磁石の平行移動は光学を変えないので、追跡が正しい) | tbende の 2 次の幾何の項。場所を特定していない |
| K1=0 の BEND の追跡のスピン (tbend + tradk) | 1 個の BEND のスピン回転角が、Emittance[]・tbendi (v2) と 2.2e-5 rad 違う (放射の分割が粗く、スライスの前のエネルギーで歳差を計算する) | K64 の精度の範囲。エネルギー損失の半分に相当 |
| B21 の K64 のソース | tcod.f の rfsw=.false. の枝は δ を動かさない (trs の 5・6 列を 0 にして 5 元で解く)。z を閉じる δ は tcod に渡す前の cod(6) で決まっているはず | temit の呼び出し側を追っていない。pysad の動作 (z を閉じる) は gs の実測と合い、物理としても妥当 |
| ANGLE=0 の BEND (tsteee・tsteer、K0 だけの補正磁石) と K0 つき MULT の端の角 | JX1 と同じ「端の蹴りを偏向角に入れるが行路は変えない」作りのはず | 端の角つきの補正磁石の例と独立な値を持っていない。直していない (ANGLE≠0 の BEND だけ直した) |
| RADINT の空洞の位相 | v2 は VOLT を単純に足す (PHI を見ない) | 位相の違う空洞のベクトル和にするかは仕様の判断 |
| tthine (長さのある六極など) の bsir0 | K64 は ∓aki·imag(cx)/al (これも次元が合わない)。pysad は実装していない (tthine は放射なしで呼ばれ、入口と出口で打ち消す) | 影響なし |

## 5. 他の担当への連絡

- `ffscmd.py::fn_track` (B33 の修正) の中の `from . import compat as _compat` が、POL + RAD (NotImplementedError で抜ける) のとき後ろの `_compat.v2()` を UnboundLocalError にしていた。
  モジュール先頭の import があるので、その 1 行を消した (k64 の `test_spin_track.py` が 5 件落ちていた)。
- カーネルが `rs[1]` (v2 の印) を読むので、`rs` は必ず `track6n.make_rs()` で作る。`dapert.py`・`solmap.py` の 1 要素の `rs` を `make_rs` に替えた (numba は範囲外を検査しない)。
- 追跡の tbendi のスピン (bsi) は B58 と同じ物理なので v2 で直した (`track6._tbendi`。track6n にスピンは無い)。

## 6. 使い方 (走る例)

```
! 矩形磁石のリングの減衰分配とエミッタンス (v2: Jx = 1.006、k64: 1.097)
MOMENTUM = 3GEV;
DRIFT D1 = (L=1.0) D2 = (L=0.5) DS = (L=0.1);
QUAD  QF = (L=0.3 K1=0.7) QD = (L=0.3 K1=-0.6);
SEXT  SF = (L=0.2 K2=0) SD = (L=0.2 K2=0);
BEND  B  = (L=1.2 ANGLE=0.2617993877991494 E1=0.5 E2=0.5);
MARK  ORIGIN = ();
CAVI RF = (L=0 VOLT=4E6 FREQ=5.0E8 PHI=0);
LINE CELLQ = (QF D2 B D2 QD DS SD D1 SF DS);
LINE RINGQ = (ORIGIN RF 24*CELLQ);
FFS;
USE RINGQ; CELL; CALCULATE; RFSW; RAD;
EMITDIV=20;                       ! 放射の分割を細かく (Jx が解析値と 1e-4 で合う)
e=Emittance[]; Print[DampingRate/.e]; Print[Emittances/.e];
RADINT;                           ! v2: 先頭が SEXT でも NaN にならない。GX・GY は dν/dδ
POL; e=Emittance[]; Print[{Polarization,PolarizationTime}/.e];   ! v2: 92.4% (平面リング)
```

`py -3.12 run_sad.py <この台本>` (v2 が既定)。K64 と同じ値は `PYSAD_COMPAT=k64`。
Python からは `emit6d.emittance6d(lat, p, rfsw=True, radcod=True, pol=True)['spin']` に `equpol` (v2 の偏極度)、`depol_rate` (1 周あたりの脱偏極率)、
`taudep` (脱偏極時間 [s])、`spin_spread` (平衡分布での n の広がり ½⟨|N u|²⟩)、`spin_orbit` (入口の ∂n/∂u、3x6)、`equpol_k64` (K64 の式の参考値) が入る。

## 7. 使える範囲と条件、使えない範囲 (v2)

| 機能 | 使える範囲 (確かめた条件) | 使えない・注意 |
|---|---|---|
| RADINT | BEND (K1 つきも)・QUAD・SEXT のリングと輸送路 (INS)。先頭が SEXT でもよい。空洞が複数あってもよい (同じ周波数)。GX・GY は一次の色収差 dν/dδ | MULT・SOL は積分に入らない (表示が出る)。BEND の K0 は入らない。結合のあるリングの V 面 (I4V) は確かめていない |
| Emittance[]・EMIT の減衰・エミッタンス | 端の角のある偏向磁石 (E1・E2)、勾配つき偏向磁石。Jx は解析値と 1e-3 (既定)、`EMITDIV=20` で 1e-4 | ANGLE=0 の BEND (補正磁石) と K0 つき MULT の端の角は K64 のまま。SOL 領域内の MULT・BEND は未対応 (元から) |
| 追跡 (RAD) の減衰 | 同上。K1 つき BEND で減衰の和が 4 | 同上。FLUC の乱数列は K64 と同じだが、端の角のある BEND では放射の大きさが変わるので値は K64 と違う |
| MomentumCompaction・BucketHeight・BunchLength (RFSW OFF) | I1/C と 2e-4 | — |
| CALC のチューン・結合 (ロール + ずれのある BEND・六極・薄い四極・空洞) | 追跡の一周写像と一致 (六極 1e-9、BEND 4e-6) | ずれた BEND の CALC と追跡の 2 次の差 (2 mm で 3e-5) は残る |
| Polarization (POL) | 平面リング 92.4%。誤差のあるリング: 脱偏極率が FLUC つきの追跡と 0.5〜1.6 倍 (3 例、§3)。`depol_rate`・`taudep`・`spin_spread`・`spin_orbit` も使える | 線形 (SLIM 相当) の見積もり。高次の側帯共鳴、運動論的偏極の項なし。SOL 領域内の MULT・BEND のスピンは未対応 (元から)。POL のときは四極を細かく分割するので Emittance[] が遅くなる |
| スピン追跡・スピン回転 (POL) の端の縦磁場 | BEND (K1 つきも)・QUAD・MULT の端。Emittance[] と追跡が 1e-5 以下で一致 | K1=0 の BEND の追跡は放射の分割の粗さで 2e-5 rad の差が残る |
| 4D の CALC の縦の Twiss 列、UnstableLevel、CalculateOptics の stabz・tracez | 上の表のとおり | — |
