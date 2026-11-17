# 47. pysad v2 の変更の一覧 (v1 = K64 に忠実な翻訳との違い)

最終更新: 2026-11-17 (v2.0.0)。

- v1 系は SAD (K64 版) のソースに忠実な翻訳で、K64 のバグや癖も gs の実動作に合わせて再現した。v2 系は、K64 のバグの修正、物理の誤りの修正、機能の追加を行う ([`../CHANGELOG.md`](../CHANGELOG.md))。
- この文書は一覧。根拠の数値・条件・使い方の詳細は作業の記録 (`history/v2_changes_*.md`) にある。記録と例の対応は文書の末尾の「記録と例」。
- K64 の側のバグの番号 (B01〜B110) は [`21_sad_bug_report.md`](21_sad_bug_report.md)。B72 以降は v2 の作業で新しく見つけたもの。
- 使った条件: 電子。3 GeV の FODO リング (24 セル 105.6 m / 91.2 m、12 セル 45.6 m / 62.4 m)、1 GeV・3 GeV の短い輸送路、0.5 GeV・2 GeV のリング、ソレノイドは 0.3〜3 T/(Bρ) 相当、
  RADLIGHT は 10.2 MeV、空間電荷は 5〜50 MeV。個別の条件は各表と記録に書いた。**実機の大きなリング・挿入光源・強い結合では確かめていない。**
- v2 の値の根拠は、解析式、独立な数値計算 (scipy の solve_ivp によるローレンツ力・Thomas-BMT の積分など)、保存則、別のコードとの一致で、gs ではない。
  gs (SAD 1.4.4.2k64、2026-10-19 ビルド) は「K64 の動作の確認」と「K64 が正しく動く機能の参照」に使った。

## (a) 互換モードの使い方

| 項目 | 内容 |
|---|---|
| モード | `v2` (既定): K64 のバグを直した動作。`k64`: v1 と同じ、K64 (gs) の動作 |
| 切り替え | 環境変数 `PYSAD_COMPAT=k64` (または `v2`)。Python からは `from pysad import compat; compat.set_mode('k64')` (`None` で環境変数に戻す)。`compat.mode()` が現在のモード、`compat.v2()`・`compat.k64()` が真偽 |
| 試験・回帰 | 既存の試験は gs との数値比較なので k64 モードで走る (`pysad/tests/conftest.py` が `PYSAD_COMPAT` の既定を `k64` にする)。v2 の試験は `test_v2_*.py` で、フィクスチャ `v2_mode` で v2 を明示する。回帰 (`run_regression.py`) も k64 で走らせる |
| gs と同じ結果が要るとき | `PYSAD_COMPAT=k64` を付ける。`py -3.12 run_sad.py lattice.sad script.sad` は既定が v2 なので、gs と数値を比べるときは必ず付ける |
| 例 | `example/v2_B/`・`example/v2_C/` は既定 (v2) で走る。K64 と同じ動作は `PYSAD_COMPAT=k64` (v2 だけの機能は k64 では動かない・無視される) |

- **モードで切り替えない変更がある**: pysad 自身の不具合の修正 (§(d))、K64 に無い機能・動かない機能のうち gs と同じ結果が出るもの (LOSSMAP、RADTAPER など)、pysad が gs に合わせて足した機能 (NOKEEPEXP、INDEX など)。表の「k64 モード」の列でモードごとの動作を書く。
- **k64 モードは「v1 と同じ」ではない点が 2 つある** (gs に合わせて変えた): (1) `For` の本体を実行しない (B71。v1 の pysad は実行していた)。(2) 利用者のコードの `BeginPackage`・`Begin` が名前空間を作る (v1 は無視して全部大域。SAD の `Packages` の読み込みは変えていない)。
- 内部の約束: 追跡のカーネルが `rs[1]` (v2 の印) を読むので、`rs` は必ず `track6n.make_rs()` で作る。

## (b) v2 で結果が変わるもの (K64 = k64 モードの動作との違い)

表の「K64」は K64 (gs) の動作で、k64 モードの動作でもある (例外は注記)。ID は [`21_sad_bug_report.md`](21_sad_bug_report.md) の B 番号。
k64 モードを持たず「従来どおり未対応」のものは「K64」の欄にそう書いた。試験名は `pysad/tests/` の下。

### 放射・エミッタンス・光学 (記録 A)

試験は `test_v2_A_*.py` (67 件)。使った条件: 3 GeV の FODO リング、偏向磁石 L=1.2 m・15 度 (ρ=4.58 m)、RF 500 MHz。

| ID | 項目 | K64 | v2 | 根拠 | 試験 |
|---|---|---|---|---|---|
| B06 | RADINT の BEND の L/ANGLE | 全 BEND に最後の BEND の L/ANGLE | その BEND 自身の L/ANGLE | 曲率の違う 2 種の BEND のリングで自前の 3x3 行列の積分と I1〜I5 が 2e-5 で一致 (K64 は I1・I4・I5 が 1e-3 以上ずれる) | `test_v2_A_radint.py::test_per_bend_radius_v2` |
| B07 | RADINT の六極 | 直前の DL を使う。L=0 は 0.2 m。先頭が SEXT だと NaN | 自分の L と K2。L=0 は寄与なし | 軌道 0 なら六極は I1〜I5 を変えない (K2=3 と 0 で 12 桁同じ) | `…::test_sextupole_first_and_thin` |
| B12 | RADINT の定数 | 単精度 (1e27、1.E-6、0.2、2.3026) | 倍精度 | ELOSS = Cγ E⁴ I2/(2π) と 1e-13 で一致 (K64 は相対 1.15e-8 小さい) | `…::test_double_precision_constants` |
| B13 | RADINT の CAVI | 最後の空洞の値 | VOLT は総和。FREQ・HARM は電圧のある最初の空洞 (違う周波数が混ざると注意) | 2.5 MV + 1.5 MV → 4000 kV。νs = √(hαV cosφ/2πE) と 2e-3 | `…::test_cavity_voltage_is_summed` |
| B72 | RADINT の `QUANTUM LIFETIME(hours)` | exp(2.3026·τ[秒]) | τ[秒]/3600 | τq = τE/2·exp(ξ)/ξ、ξ = (δmax/σδ)²/2 | `…::test_quantum_lifetime_hours` |
| B73 | RADINT の `LINEAR CHROMATICITY` | −I9/4π (磁石の中だけの ∮γ ds で色収差にならない。例: −2.36、実際は −7.30) | 運動量を ±1e-4 ずらした光学のチューンの差分 dν/dδ (六極・偏向の寄与を含む) | 同じスクリプトの `DP0` をずらした CALC の差分と 2e-5 で一致 | `…::test_chromaticity_is_tune_derivative` |
| B74 | RADINT の MULT・SOL | 黙って無視する | 場のある MULT・SOL があれば「入っていない」と表示 (積分には入れない)。k64 モードは表示なし | — | `…::test_mult_warning` |
| B75 | Jx の端の項 (docs/46 の「RADINT と EMIT の Jx の差」) | Emittance[]・EMIT・追跡は、端の角 e の蹴りを偏向角 θ に入れるが行路 al1 を変えない。I4 の端の項が 2 倍になる。矩形磁石のリングで Jx = 1.097 | 端の角のある面で場の中の行路 al1 を x·tan(e) だけ短くする (微分も D 行列に入れる) | (a) 解析式 Jx = 1 − I4/I2 と (b) 厳密な幾何 + Robinson の Je = 2 + dln l/dδ が 2e-6 で一致 (E1=E2=0: 0.92082、0.1309: 0.94147、0.5: 1.00575)。**RADINT が正しく K64 の EMIT が誤り**。v2 の Emittance[] は 1.00586 (EMITDIV=20 で 1.00580)、追跡の一周写像の固有値は 1.0033。gs の Emittance[] の Jx は 1.0974、同じ gs の RADINT は 1.00575 | `test_v2_A_damping.py` (7 件) |
| B76 | K1 つき BEND の追跡の曲率半径 | tbendi は放射で運動量が変わっても rhoe を入口の値のまま。減衰の和が 4 にならない (K1=0.1 で 4.05、Jx=0.7385) | スライスごとの放射のあとに rhoe = rhob·(1+δ) を更新 | 和が 4.000、Jx = 0.6875 (解析 0.68762、RADINT 0.68762、Emittance[] 0.6883) | `…::test_track_edge_and_sum_rule_v2` |
| B77 | MomentumCompaction (Emittance[]) | RFSW OFF は −R56/C だけ (Σ R5j η_j が抜け、FODO の例で 0.022898、正しくは 0.022137: 3.4% 大きい)。RFSW ON は μ² から (μ²/12 だけ大きい: 0.25%) | RFSW OFF は −(R56 + Σ R5j η_j)/C。RFSW ON は 4 sin²(μ/2) | 自前の I1/C と 2e-4 (FODO の例は 6 桁)。BucketHeight 0.009877 → 0.010045 (RADINT 0.010045、RFSW ON 0.010047)。docs/46 の「BUCKET HIGHT の差」はこの α の差 | `test_v2_A_radint.py::test_momentum_compaction_is_i1_over_c_v2` |
| B09 | tchge の入口のロール + ずれ | tbrote が回転行列を上書きし、一周行列に回転が入らない | 入口の回転を行列に掛ける | `_tchge_jac` が `_tchge_cod` の数値微分と 1e-8 で一致。ロール + ずれのある六極・BEND のリングの CALC のチューンが追跡の一周写像の固有値と一致 (六極 1e-9、BEND 4e-6。K64 は 0.02〜0.18 ずれる) | `test_v2_A_tchge.py` (5 件) |
| B56 | tbrote のスピン回転 srot | 行の逐次更新で行列積にならない | (I+rr)·srot | 直交性が 1e-14 で保たれ、入口と出口が互いに逆 (K64 は 1e-5 以上崩れる) | `test_v2_A_tchge.py` (3 件) |
| B11 | 4D の CALC の縦方向の Twiss 列 | 作り直さず 6D の値が残る | 既定値 (AZ=0、BZ=1、他 0) | 4D では縦方向を計算しない | `test_v2_A_ffs.py::test_calc4d_resets_longitudinal_columns_v2` |
| B15 | `UnstableLevel` | 常に 0 | 不安定な点の数 nstab | MatchingResidual = 10·nstab + NetResidual (例: 2 と 20) | `…::test_unstable_level_v2` |
| B16 | CalculateOptics の stabz・tracez | 未初期化 (k64 モードは 1 と 0) | idim=3 は tracez = 2cos(μz) と stabz、idim=2 は stabz=1・tracez=0 | Emittance[] のシンクロトロンチューンと 2e-3 | `…::test_calculate_optics_idim3_returns_z_trace_v2` |
| B39 | tradke (SOL 内) の dbeam(8) | B12 − B44 + b·B13 | B12 − B34 + b·B13 | 力学的運動量 (px + b y, py − b x) の共分散 A B Aᵀ | `test_v2_A_spin.py::test_kinetic_momentum_covariance_formula`、`…::test_tradke_px_py_diffusion_uses_b_y_py_v2` |
| B58 | BEND の端の縦磁場の積分 bsir0 (Emittance[]・追跡の POL) | 出口の項が常に捨てられる。K1 の項の次元が合わない。スピン軌道結合に端の項の y・x 依存が入らない。追跡の tbendi は y/ρ がない | 入口で +(y/ρ + k x y)、出口で −(…) (k = K1/L)。軌道についての微分もスピン軌道結合に入れる。追跡の tbendi も同じ | ∫Bs ds/(Bρ) = ψ(x,y)。軌道のずれた 1 要素のスピン回転が Emittance[] と追跡で一致: K1=0.3 の BEND 1.7e-6 (K64 2.0e-3)、K1=0 の BEND 7.5e-5 (K64 1.0e-3)。スピン軌道結合は追跡の差分と 2% (K64 は 1/6) | `test_v2_A_spin.py` (3 件) |
| B64 | MULT の bsir0 (Emittance[] の POL) | imag(cx)/aln (次元 1/m)。四極は入れない。分割は放射の大きさだけで決まる | ψ = Im Σ k_n z^(n+1)/(n+1)! を入口で足し出口で引く。四極は k x y。POL のときは四極・MULT を k·(L/n)² ≤ 2e-3 になる数に分割 | MULT (K2=30) 3.1e-7 (K64 9.4e-5)、MULT (K1・K2・SK2) 3.8e-6 (3.1e-4)。四極の分割が 1 個だと歳差の微分が追跡の差分と 50% 違い、細かくすると 2% で合う | 同上 |
| B55 | `Polarization` (POL の平衡偏極度) | NaN (spintrm の exp(cb²−cφ)·cerfc(cb) が 0·∞)。k64 モードも NaN | (1) exp(cb²)·erfc(cb) を Faddeeva 関数 w(i·cb) で。(2) 平衡偏極度は線形のスピン軌道結合 (DK の (11/18)\|γ∂n/∂γ\|² の項) の脱偏極率 rate から P = (8/5√3)·sdamp/(sdamp+rate) | 平面リングで 92.376% (Sokolov–Ternov)。誤差のあるリングは FLUC つきのスピン追跡で測った脱偏極率と 0.5〜1.6 倍で合う (3 例: 3 GeV 1.6 倍、2.903 GeV 1.2 倍、4.2 GeV 0.53 倍。16 粒子、統計 ±25%)。K64 の式 (桁あふれだけ直す) は 3 桁以上大きい | `test_v2_A_spin.py` (5 件) |

### 追跡・要素 (記録 A・C・C2・C_measure)

試験は `test_v2_A_track.py`、`test_v2_C_*.py`、`test_v2_C2_*.py`。使った条件: 電子、0.5・1・3 GeV、周長 91 m の FODO リングと 4〜6 m の輸送路、ソレノイド 1〜3 T (bz = 0.3〜3 1/m)。

| ID | 項目 | K64 | v2 | 根拠 | 試験 |
|---|---|---|---|---|---|
| B08 | tsolqur (RAD、SOL 領域内の四極) | 最後の半ステップで最後の粒子の bzp を全粒子に使う | 粒子ごとの BZ/(1+δ) | 1 粒子の結果が、一緒に追う粒子に依らない (差 0。K64 は 2.6e-6) | `test_v2_A_track.py::test_solenoid_quad_rad_is_independent_of_other_particles_v2` |
| B14 | tsolqum (放射なしの MULT) の z の dv の項 | −0.25·(ap−dv)·aln0 | tsolqu と同じ −dv·aln0 | 同じ K1・L の QUAD と MULT の z が 1e-12 で一致 (K64 は 20 MeV・L=0.5 m で 1.9e-5 m ずれる) | `…::test_mult_z_equals_quad_and_drift_v2` |
| B10 | tdusrn の NaN の anp | npmax = 10000 回まわり乱数を 20001 個引く | NaN・±Inf なら光子を出さず乱数は 1 個 (有限の anp は K64 と同じ乱数列) | 同じ種で同じ戻り値と種 | `test_v2_C_track.py::test_b10_*` (3) |
| B17 | NaN と MAX/MIN | gfortran の MAX/MIN が NaN を無視する前提 | 乱数の消費だけ B10 で NaN を明示的に扱う。asinz・gfloor・sqrt1 は変えていない | 要判断 1 (§(e)) | (B10 の試験) |
| B33 | セッション最初の TrackParticles の最初の tapert | l_track = 0 で、LOSSAMPL の外の粒子が落ちず追跡が続く | 最初の tapert で落ちる (座標は 3·LOSSAMPL に頭打ち、落ちた周 1、落ちた要素は始点の 1 つ前) | 2 回目以降の K64 の動作と同じ | `test_b33_*` (2) |
| B35 | SOL 領域内に DRIFT BEND QUAD SOL MULT CAVI MAP APERT MARK MON 以外の要素 | SIGSEGV (k64 モードは NotImplementedError。L=0 は素通し) | `SolRegionError` (K64 の tsol.f が出すはずの文)。長さ 0 の SEXT なども黙って通さない | 仕様 (tsol.f の case default) | `test_b35_*` (1) |
| B37・B78 | UND の磁場 | 磁場のキーワードが無く、ドリフト。さらに undulator.f の陰的刻みの右辺で dpds を二重に掛け (B78)、磁場があっても収束力が刻み幅に比例して消える。k64 モードは従来どおり (FBX は NotImplementedError、光学は恒等) | `FBX`・`FBY` [T] を受ける (追跡・CALC・Emittance[]・EMIT。RAD で放射も)。磁場つきカーネルは dpds を 1 回だけ掛ける。**v2 では磁場なしでも UND は長さ L のドリフトとして光学に入る** | 平面アンジュレータの解析式: 軸上の軌道 (出口で x = x' = 0)、行路長 (1e-6)、垂直の収束 k = F0²/2 (5e-4)、シンプレクティック性 1e-7、放射損失 (2/3) r_e γ³ (F0²/2) L (1e-4)。リング: Δνy = ∫βk ds/4π (1%)、ΔU0 (2e-3)、減衰率の和 = 2U0/E (1e-3) | `test_v2_C_und.py` (9) |
| B38・B79 | PHSROT | 行列がすべて 0 で粒子の座標が 0 になる。PX PY PZ は指定できない。拡散の乱数が非対角の U で共分散を壊す (B79)。k64 モードは追跡が NotImplementedError、光学は恒等 | 意図された写像 M = N·Rot(2πPX, 2πPY, 2πPZ)·(1−D)·N⁻¹ + 拡散 (RAD・FLUC)。拡散は N U (√λ ⊙ ξ) | M0 のシンプレクティック性 1e-13、固有値 = exp(±2πi P)、リング (ST PHSROT) の Twiss・チューン = 入力 (1e-9)、平衡エミッタンス b/(2d−d²) (Emittance[] 1e-6、追跡の統計 8%)、減衰率 = D | `test_v2_C_phsrot.py` (5) |
| B40・B62・B99 | DELCOR | MARK には DETR・DPX、CAVI には HARM・PHI に (符号も逆に) 書く (B40)。corlength≠0 は波数 −1 の成分 cv(n) が未設定で、負の波数の重みが 1 つずれる (B62)。corlength≠0 の規格化が複素 c の絶対値の 2 乗平均 = 1 で DX の rms が delx/√2 (B99) | MARK・CAVI も DX・DY に他の磁石と同じ符号で入れる。cv(n) を対称に置く (乱数の順は同じ)。`s = 1/sqrt(2Σw²+1)` で DX の rms = delx、DY の rms = dely | 相関長 ≫ 周長で全要素のずれが同じ (QF・RFC・ST が 5% 以内、5 km で全桁同じ)。相関関数 exp(−Δ²/(2 x0²)) (300 標本で 5%、虚部 2%)。場の rms (x0 = 15 m): v2 0.987・1.009、k64 0.708・0.703。DELCOR コマンド (62.4 m のリング、150 回): rms(DX)/delx = 1.021 (k64 0.713)。B99 は作者の意図の裏づけが無い (推定) | `test_v2_C_delcor.py`、`test_v2_C2_delcor.py` (7) |
| B65 | RADLIGHT の軌道 | 横方向の座標が入らず、BEND で向きが変わらず、BEND の中を記録しない | 横の軸を幾何から取り、BEND は円弧で進め、中も 1/γ の角度あたり `RADLIGHTDIV` 点 (既定 16) 記録。曲げのある MULT は出口だけ、CAVI も記録 | 振幅 0 の粒子の点 = 幾何 `LINE["GX"/"GY"/"GZ"]` (1e-11)、時刻 = 行路長/(βc) (1e-9)。偏向磁石の放射場のピーク E = e β²/(4πε0 ρ R (1−β)²) と RadiationField が 0.25% (RADLIGHTDIV 16)、5e-5 (64) | `test_v2_C_radlight.py` (3) |
| B69・B80 | SOL 領域内の BEND・MULT・CAVI・APERT | gs が領域内の BEND の偏向をほぼ無視する。k64 モードは NotImplementedError。放射 (tdrift) は 1 刻みの弦から損失を出すので、らせんの回転角 bz·L が大きいと小さく出る (相対 (bz·L)²/12。bz·L = 1 で 8%) (B80) | BEND = 直線の座標のままの一様な双極磁場 (tsol.f の `tdrift(al, bzs, phiy, phix)`)、MULT = tmulti1、CAVI = tmultiacc を bz つきで。領域内の BEND・MULT は bz·L/ndiv ≤ 0.02 に刻む。APERT で落ちる | ローレンツ力 dp/ds = (p × B)/pz の数値積分 (solve_ivp、rtol 1e-13): BEND・K0 だけの MULT 1e-14、K1 つき 1e-10、六極・八極・スキューは 1e-8。放射損失 (2/3) r_e γ³ (K0/L)² L が bz = 0・0.4・2 で 2e-4。bz → 0 の領域内の要素 = 領域の外 (リングの Twiss・Emittance[] が 2e-7)。CALC の軌道 = 追跡 (1e-12)。bz = 3 T のリングで減衰率の和 = 2U0/E (2e-3) | `test_v2_C_sol.py` (30) |
| B81 (N6) | スピン追跡 sprot の縦磁場 bsi | bz = bz0/pzm + bsi (p0 で規格化したまま z 軸の向きに足す)。δ≠0 で (1+a)B∥ の歳差が 1/(1+δ) 倍にならない (δ = −5e-3、0.12 rad で 6e-4 rad)。傾いた軌道で bsi·ẑ の垂直成分が B⊥ に二重に入る。双極 + ソレノイドの bsi が T-BMT と 3.6e-3 rad 合わない | SOL 領域と UND では ∫B∥dl/(Bρ) = [(bz pz + Bx px + By py)·al/pz + ψ]/(1+δ) を運動量の向きに置く | Thomas-BMT の積分: v2 1e-7〜1e-5、K64 の式 6e-4〜4e-3 | `test_v2_C2_spin.py` (17 件中) |
| B82 (N7) | tsolqur の縦磁場と放射の基準 (SOL 領域内の四極、RAD) | bsi を最初・最後のスライスは半分で上書き (1 スライス分足りない)。放射の基準 pxr0 が正準運動量のままで、最初のスライスに (bz/2)(y, −x) の見かけの偏向 | スライスごとに bz·aln/(1+δ)。放射の基準は力学的運動量 (高速追跡も) | 軸上の粒子のスピンの回転 (1+a) bz L: K64 は半分、v2 は 5e-15。軸外の放射損失: v2 の四極と MULT が 3e-4 で一致、K64 は 2.2 倍 | 同上 (`[quad_sol]`) |
| B83 (N8) | sprot の刻みの歳差の角 (すべての要素の追跡のスピン) | (1+aγ)·2 sin(θ/2)。刻みが粗いと θ²/24 小さい (0.26 rad の BEND を 3 刻みで 3e-4) | 厳密な θ | 平面リングの 1 周の回転角 = 2π aγ (K64 の式は 4e-3 rad ずれる) | `test_ring_spin_tune_with_solenoid[0.0]` |
| B85 (N11) | AUTOPHI | TRPT でないとき phis = phic で dphis = 0 になり、DZ·w も PHI も効かない。光学 (phii = phic) と合わない。全空洞が AUTOPHI だと実効の RF の傾きが 0。k64 モードは NotImplementedError | 「位相の基準 = 設計軌道の到着時刻」(TRPT のときの形 phis = 0)。追跡 (領域の内外)・CALC (TRPT)・Emittance[] に実装。傾きは軌道の位相での値を足す | 設計軌道が DZ = 3 mm 早く着く輸送路: ありは設計軌道の粒子の δ が 0 (CALC 1e-12、追跡 1e-9)。gs (TRPT) と CALC の DZ・DDP・追跡が全桁一致 | `test_v2_C2_solrf.py` (AUTOPHI の 2、Emittance[] の 1) |
| B86 (N13) | SOL の trots の z の DZ | 入口は z に DZ だけ、出口は DZ·(1 − dv) | 入口・出口とも (1 − dv)。k64 モードは SOL のずれ・回転が未対応 | 差は dv·DZ (3 GeV で 1e-13 m 以下) | — |
| B41・B87・B88 (MEASURE の初期分布) | MEASURE の初期分布 (ttinit) | 軌道の揺らぎを MARK の添字 26..29 (ZPX ZY ZPY DP) から読む (B41)。QUIET の重心合わせの基準が分布の作り方と合わない (B87)。TWAKE/LWAKE ON で全粒子 z = 0 (B88)。MEASURE は無出力 (B42)。k64 モードは MEASURE を無視 | JDX JDPX JDY JDPY (34..37) を読む。QUIET は分布と同じ写像で求めた重心に合わせる。z の分布を作る (`ttinit(b41=True, center_fix=False, wake=True)` で K64 の動作) | n=200000 の JITTER ON の平均が v2 の基準と 4σ 以内 (K64 の基準から 20σ 以上)。MK の例で px の差 2.5e-5 (統計誤差 3e-7)。B41 は MARK のキーワード表の添字。**ソースの読みのみ (gs は MEASURE が無出力)** | `test_v2_C_measure.py::test_b41_*`、`…::test_ttinit_quiet_centroid_*` |
| B42 (MEASURE・TRCOD) | MEASURE・TRCOD・LOGMEASURE・MAIN の `TRACK TURNS=-1` | MEASURE は無出力、TRCOD は SIGSEGV、MAIN の TRACK は落ちる。k64 モードは MEASURE・TRCOD を「未対応 (無視)」、TURNS=-1 は 1 周のリングの追跡 | K64 のソースの意図どおりに動く (§(c))。MEASURE の間だけ TRPT を ON にして戻す | 初期分布の 2 次モーメント (n=100000) が統計誤差 1.3% 以内、tbdecoup の非結合 emittance が 4×4 行列の固有値と 1e-10、MEASURE の出口の分散行列 M Σ0 Mᵀ と 1.9e-7、TRCOD の軌道が CALC (INS) と 5.4e-11。**gs との比較は無い** | `test_v2_C_measure.py` (19) |
| B100〜B105 | TRCOD の FIXCOD・MVCOD・STNCOD・ALIGN、ORBITJITTER、PUTMEASURE・GETMEASURE | ソースの誤り・旧い番号の名残 (docs/21 の B100〜B105)。K64 の FFS からは基準を読めない。k64 モードは語として認識しない / 未対応 | §(c) の機能として実装 (基準は REFERENCE と同じ参照光学に置き `Twiss["RDX"]`・`Twiss["DDX"]` で読む) | FIXCOD の差 Δx が解析式と相対 1e-5 (輸送路)・1e-4 (リング)、ALIGN が平行移動と 1e-13 (BEND を含むと 4e-10)、PUTMEASURE → GETMEASURE が 7 桁。**gs との比較は無い** | `test_v2_C2_measure.py` (20) |

### 集団効果 (記録 C・C2・C_spch・C2_spch)

| ID | 項目 | K64 | v2 | 根拠 | 試験 |
|---|---|---|---|---|---|
| B34 | BEAMBEAM の傾き回転 | cost = sin(傾き角)、sint = 傾き角 (回転にならず座標が潰れる) | (cos t, sin t) の回転。強ビームに x-y 結合があるときだけ結果が変わる | 回転に対する共変性: t = 0.3 で 1e-9 (蹴りに対する相対) | `test_v2_C_beambeam.py` (2) |
| B84 (N10) | 空洞の wake の重み (tcav.f・tmulti.f の `fw=fw0*(wsn+ws(n-1)*.5d0)`) | 空洞が 2 つ以上に分割されると蹴りの合計が (1.5n − 0.5)/n 倍。k64 モードは再現 | 台形の重み fw0·(ws(n) + ws(n−1))/2 (合計 1)。領域の外の空洞・RF つき MULT も | 位置によらない表で、エネルギーの変化が分割数によらない: 20 MV (1 分割) と 200 MV (4 分割以上) で 5e-3。K64 は 1.25〜1.5 倍 | `test_v2_C2_solrf.py::test_wake_weight_sum_is_one_in_v2_and_one_and_half_in_k64` |
| B66 | 旧式の `WAKE` コマンド | 配布に無い WAKECOMMAND.n を読もうとして失敗 (tfwake.f は死んだコード) | `WAKE L 要素 s1 w1 … ENDWAKE` (横は `T`) を WakeFunction の表として登録 (CAVI だけ。行をまたいで書ける §(c)) | 同じ表を `WakeFunction[…] = …` で与えた追跡と 2e-9 (gs の参照値と一致済みの追跡) | `test_v2_C_wakecmd.py` (3)、`test_v2_C2_wakecmd.py` (6) |
| B59 | CSRSetupOY の戻り | dpjn が zjn と同じ記憶域 (zj に運動量が入り、dpj は 0)。点が無いと i1 が未初期化 | zj は位置、dpj は運動量。点が無ければ失敗 (CSR::irreg) | 蛙跳びの関係 z_{k+1} − z_k = p_k r56 dt1 (1e-9) | `test_v2_C_csr.py` (2) |
| B89〜B92 (SPC1〜SPC4) | SPAC (DRIFT の `RADIUS>0`) の空間電荷 | 力が 1/γ 倍に弱い (B89)。px を二重に (1+δ) で割る (B90)。メッシュが整数の割り算で np<1000 では縦 10 枚、np=2000 でも横 8×8 (B92)。k64 モードは SPAC・PSPAC・WSPAC を追跡に入れない | バンチの静止系の 3 次元のポアソン方程式 (円筒の壁は像の電荷) で正しい係数。24〜64 ノードの 2 次 B スプライン | 一様な長い円柱 x'' = K x/a² (K = 2 r_c λ/(β²γ³))、KV 包絡方程式 (5 MeV・1e11 個で a(1 m) = 3.996 mm vs 追跡 2·rms 5e-3)、Bassetti-Erskine 5e-3、リングのチューンシフト (Laslett) 追跡/モデル 0.97〜1.00。gs は 5 MeV で rms 0.6226 mm (解析 K: 1.639、K/γ: 0.645、v2 1.623) | `test_v2_C_spch.py` (22)、`test_v2_C2_spch.py` (35) |
| B91 (SPC3) | MULT の `RADIUS` の空間電荷 | spkick の呼び出しの引数が合っていない。最後の回だけ長さが半分 | ANGLE=0・VOLT=0・SOL 領域の外の MULT に入れる。四極を MULT で書くと全周に入る | 四極の中の KV ビームの包絡 2·rms 0.07%・0.2%、リングのチューンシフト 追跡/モデル 0.94〜0.97 (SPAC)・Laslett の式と 1% | `test_v2_C2_spch.py` |
| B93〜B98 (SPC5〜SPC10) | PSPAC (SPCH 要素) と SPAC の小粒子数 | `PSPACNZ=1` (既定) で gs が SIGSEGV (B93)。λ(z) の配分の不連続 (B94)。周の終わりのキックの長さ (B95)。状態が static (B96)。tspac は死んだコード (B97)。SPAC の np が少ないと SIGSEGV・Bus error・NaN (B98) | `PSPACNZ<4` はバンチに合わせた 64 ノードを自動で張る。2 次 B スプライン。区間の始めから数え終わりに残りを入れる。呼び出しごとに作り直す。粒子数によらず走る (2・5・15・33・81 個) | PSPAC の 3 周後の rms が gs と 4 桁目まで (1.777/1.042 mm と 1.776/1.042 mm)、Laslett の式 (PSPAC: 追跡/モデル 0.98・1.004)。B98 は gs で np = 2〜300 を走らせて原因を特定 | 同上 |

### 特殊関数・数値 (記録 B。pysad 自身の不具合も含む)

試験は `test_v2_B_specfn.py` (13 件)。実数の倍精度。mpmath は無く、scipy.special と scipy に依らない独立な計算 (閉じた式、級数、AGM、積分表示、Wronskian) の 2 つを突き合わせた。

| ID | 項目 | K64 | v2 | 根拠 | 試験 |
|---|---|---|---|---|---|
| B43 | `Hypergeometric2F1`・`BesselJ` | `Hypergeometric2F1[.5,.5,1,.99]` = 2.348414、`BesselJ[.5,30.]` = −.143924709693 (打ち切り) | 2.352715816780、−.143929653370 (pysad は以前から真値。モードによらない) | 2F1(½,½;1;m) = 2K(m)/π を AGM で、J_{1/2}(x) = sqrt(2/πx) sin x。相対 1e-13 | `test_b43_hyp2f1_besselj_precision` |
| B44 | `ExpIntegralEi` | `Ei[1.5]` = .164−2.985 I、`Ei[-.1]` = −2.2265 | 3.301285449130、−1.822923958419 (実数。モードによらない) | 級数 γ + ln\|x\| + Σ x^k/(k k!)、負の側は −E1(x)、大きい x は漸近級数。1e-10〜2e-12 | `test_b44_b54_expintegralei` |
| B45 | `HermiteH[3,.7]`・`LegendreP`・`Binomial[-3,2]` | NaN | −5.656、+.274954541697 (`LegendreP[3,1,.4]`)、2.73 (`LegendreP[2,2,.3]`)、6 | H3 = 8x³−12x、P_3^1 = −(1−z²)^½ (15z²−3)/2、P_2^2 = 3(1−z²)、P_4^3 = −105 z (1−z²)^{3/2} | `test_b45_*` |
| B46 | `HypergeometricU`・`BesselY`・`BesselK`・`BesselJ[-1,…]`・`BesselJZero` | INF・NaN・無応答 | 2、−1.534238651350、1.806573612779、−.557936507910、2.404825557696 (モードによらない) | U(1,2,z) = 1/z、Y_0 の級数、K_{3/2} の閉じた式、Wronskian J_1 Y_0 − J_0 Y_1 = 2/(πx) | `test_b46_*` |
| B47 | `InverseErf` | `InverseErf[.3]` = .272462711811 (Newton の係数が sqrt(π/2)、相対 1e-8 の誤差) | .2724627147267543 (scipy の erfinv)。`InverseErfc` は小さい z で桁落ちしない erfcinv。**モードで変わる** (k64 は K64 の値) | erf(InverseErf(x)) = x が 4e-16 以内 | `test_b47_*`、`test_fitgaussian_uses_true_inverseerf_in_v2` |
| B53 | 大きい引数の `BesselJ/K/Y`・`FresnelS`・`AiryAi` | `BesselJ[0,33.]` 相対 3e-3、`FresnelS[6.1]` = 1.3e7、`AiryAi[15.5]` = −16 | 真値 (モードによらない) | J_n・I_0 の積分表示、K_ν の積分、Fresnel の Simpson、Ai = (1/π) sqrt(x/3) K_{1/3}(2/3 x^{3/2})、Ai Bi' − Ai' Bi = 1/π | `test_b53_*` |
| B54 | `PolyLog` | 複素数 (虚部 1e-15 の雑音つき) | 実数 (z ≤ 1)。z > 1 は切断線の下側 (Li_2(2) = π²/4 − iπ ln 2)。非整数次は \|z\| の上限なし (反転公式・積分表示。s ≤ 7.9、\|z\| ≤ 1e12 で相対 1e-11) | Li_2(½) = π²/12 − ln²2/2、Fermi–Dirac 積分、Li_{s−1} = z dLi_s/dz | `test_b54_polylog_*` |
| B18 | `PolynomialFit` (x の平均が 0) | `xl=1` のまま二項係数の和を取り係数が壊れる (`{{-1,1},{0,0},{1,1},{2,4.1},{-2,4}}` の 2 次で {1.0257, 2.0486, 1.0143}) | {−.0085714, .02, 1.0142857}。**モードで変わる** (k64 は K64 の値。平均が 0 でないデータは両モード同じ) | numpy.polyfit と 1e-13、3 次の多項式の係数の回復 | `test_b18_polynomialfit_mean_zero_v2` |
| B52 | `Fit[data,fun,{x1,x2},…]` | narg を出して最小にならない値 | 最小二乗解 (モードによらない) | numpy.linalg.lstsq と一致 (係数 1e-6、ChiSquare 1e-8) | `test_b52_*` |
| B61 | 複素の初期値の `FindRoot` | 壊れた値 (`{z->Restrict[(.25+I),…]}`)。k64 モードは TypeError | 複素の Newton 法 (§(c)) | z² = −1 → i、e^z = 2 → ln 2、z³ = 1、連立 zw = i・z + w = 2 の残差 1e-12 以下 | `test_b61_findroot_complex_v2` |

### スクリプト・FFS (記録 B)

試験は `test_v2_B_script.py` (25 件、うち 1 件は POSIX だけ)、`test_v2_B_ffs.py` (14 件)。使った条件: 1 GeV・周長 6〜7 m の四極だけのリング、3 GeV・周長 105.6 m の FODO 24 セル。

| ID | 項目 | K64 | v2 | 根拠 | 試験 |
|---|---|---|---|---|---|
| B71 | `For[ini,cond,incr,body]` | 本体を実行せず incr だけ繰り返す。3 引数は未評価。**k64 モードは gs と同じにした** (v1 の pysad は本体を実行していた) | 本体を実行 (Break・Continue も効く)。3 引数は本体なしのループ | t = Σ i² = 30 | `test_b71_for_body_v2`、`test_b71_for_body_k64` |
| B70 | `MessageName[f,"usage"]` (関数形) | f を無視して `Undefined MessageName General::"usage"` | f の本文を返す (無ければ `General::tag` の本文)。モードによらない | `General::mytag` の継承は gs の `h::mytag` の実測と同じ | `test_b70_messagename` |
| B22 | `MON`・`MCOD`・`TCOD` | コマンドでなく、後ろの語が式として読まれて `Missing value for PM*` で文が止まる | 「K64 に無いコマンド」と警告して文全体を読み飛ばす (定義はしない。k64 モードも同じ) | K64 の tffsa.f にこれらのコマンドは無い | `test_b22_mon_mcod_tcod` |
| B24 | `LINE["LENGTH","$$$"]` | narg エラー | エラーのまま。メッセージに位置の取り方 (`LINE["POSITION",loc]`、`LINE["LENG",loc]`・`LINE["S",loc]`) を足した | K64 のソース (tftwiss.f:508) が `LENGTH` を引数 1 個に限っている (仕様) | `test_b24_line_length_with_location` |
| B26 | `ORG B3 1.5 -2.5 0.5 10 20 30` | `1.5-2.5` と読まれ値が 1 つ足りない | 空白で区切った `-2.5` は別の値。`3/2`、`2 *5` などの式は書ける | K64 で曖昧さのない書き方 (`1.5 0-2.5 …`) と同じ ORG の値 | `test_b26_org_negative_numbers` |
| B30 | MAIN の `OFF RFSW` | FFS に届かない (FFS は RFSW ON で始まる) | 届く (`ON RFSW` も)。FFS の中の `RFSW;`・`NORFSW;` が後から優先 | FFS の `NORFSW;` と全桁同じ追跡結果 | `test_b30_main_off_rfsw_reaches_ffs` |
| B67 | MOMENTUM の既定 | 30 GeV で、書き忘れても何も出ない | 既定は同じ 30 GeV。MOMENTUM を MAIN でも FFS でも設定しないまま運動量を使う計算 (追跡・EMIT・Emittance・RADINT・CALC6D・DA ほか) に入ると 1 度だけ警告 | 明示の `MOMENTUM=30e9` では警告しない | `test_b67_momentum_warning` |

## (c) 新しく使える機能と、モードによらない追加・修正

「モード」の列: **v2** = v2 モードだけ (k64 モードは従来どおり未対応・無視・恒等)、**両** = 両モードで同じ (pysad の追加・修正、または gs と同じ結果が出る機能)。
例は `example/v2_B/`・`example/v2_C/` (README に表がある)。使える範囲は「確かめた条件」で、それ以外は未確認。

### 放射・光学・要素

| 機能 | モード | 使える範囲と条件 | 使えない範囲 | 例・記録 |
|---|---|---|---|---|
| **UND の磁場** (`FBX`・`FBY` [T]。B37) | v2 | 追跡 (`TrackParticles`、MAIN の TRACK、DA)、CALC の Twiss、Emittance[]・EMIT・CALC6D・Touschek。RAD で放射 (FLUC も)。`SLICE` は `POLE` の 20 倍以上 (刻みが周期の整数倍だと磁場が見えない)。`KX`・`QY` は K64 と同じ。試験は 3 GeV、1 T・1.8 T、λ = 5 cm、20 周期。磁場なしでも UND は長さ L のドリフトとして光学に入る。`SLICE=0` は恒等 | `BX`・`BY` という名前 (K64 にも無い。ValueError)、端の磁極の整形。PHI≠0 の入口はベクトルポテンシャルの飛びが蹴りになる (K64 の式のまま) | `example/v2_C/C1_undulator.sad`、`history/v2_changes_C.md` §2.1 |
| **PHSROT** (B38) | v2 | PX PY PZ はチューンの単位の位相、Dij は正規座標での減衰行列 (RAD で効く)、Bij は 1 回あたりの拡散行列 (RAD かつ FLUC)。R1..R4、EX EPX EY EPY ZX ZPX ZY ZPY も使える。追跡 (track6n)、CALC、Emittance[]・EMIT・CALC6D・Touschek。(z, δ) の回転は向きが逆 (matrix.f の m_prot_set_6d。縦のチューンは −PZ と出る)。分散つきの Twiss を CALC で読むと BX・EX は入力と 1e-3 ほど違う | k64 モード。拡散行列に負の固有値があると ValueError。スピン追跡の経路 (track6) では位相空間の写像だけで、PHSROT はスピンを動かさない (リング全体を 1 個で表すとスピンチューンは 0。偏極の計算には使えない) | `C2_phsrot.sad`、C §2.2、C2 §2.3 |
| **SOL 領域内の BEND・MULT・CAVI・APERT** (B69) | v2 | BEND: `L`、`ANGLE + K0`、`ROTATE + DROTATE`、`RAD`。直線の座標のままの一様な双極磁場。MULT: K0..Kn と SK0..SKn、`DX DY ROTATE`、端のフリンジ、RAD、`VOLT`。CAVI: `L`、`VOLT`、`PHI`、`DX DY ROTATE`、端のフリンジ。追跡、CALC (INS/CELL) の Twiss・軌道、Emittance[]・EMIT。条件: bz = 0.3〜3 1/m で確かめた | BEND の `K1`・`DX`・`DY` (NotImplementedError。四極成分は QUAD か MULT で)、ANGLE≠0 の MULT、`V1 V20 V11 V02`、SEXT・COORD・BEAMBEAM など (SolRegionError)、CALC6D・RADCOD つきの CALC (SADError)。未確認: 領域内の MULT・CAVI の POL、RF つき MULT の RAD、bz ≥ 3 1/m を超える場 | `C3_sol_region.sad`、C §2.3 |
| **ずれて傾いたソレノイド** (`GEO=1` の SOL の `DX DY DPX DPY`、`DX DY DZ CHI1..3`、端の放射 `F1`) | v2 | 追跡 (高速追跡と track6)、CALC の Twiss と軌道、Emittance[]・EMIT、幾何 (`LINE["GX"…]`) と幾何の fit。領域の中の要素の軸はソレノイドの軸 (領域の中の Twiss・軌道はソレノイドの系)。`F1≠0` は RAD で端の蹴りを長さ F1 の一様な偏向として放射させる。領域内の BEND・MULT・CAVI・K1=0 の四極も軌道を追う | CALC6D・RADCOD つきの CALC、追跡を領域の途中から始める・終える、領域より上流に軸まわりの回転 (GCHI3≠0) があるラティスで追跡を先頭以外から始めた場合、領域内の APERT のスピン追跡。未確認: 軸が設計軌道と 10 mrad 以上ずれる場合、GEO のない領域に手で書いた CHI の符号の約束。gs との数値比較は行っていない | `C2a_sol_tilt.sad`、C2 §2.1 |
| **領域内の空洞の TRPT の加速・wake・AUTOPHI** | v2 | `TRPT; RFSW;` で領域内の CAVI・`VOLT≠0` の MULT が加速する (追跡、CALC (INS))。`LINE["GAMMABETA", k]` が各位置の基準運動量。wake は領域内の CAVI・RF つき MULT に付けられる (TRPT のとき)。`AUTOPHI=1` は空洞の位相を設計軌道が着く時刻 (Twiss の DZ) を基準に測る。RF つき MULT の Emittance[] は bz≠0 の領域の中だけ | スピン追跡の経路 (track6) での TRPT の加速、`V1 V20 V11 V02`、ANGLE≠0 の MULT。AUTOPHI の意味は K64 のソースの読み (SADHelp に説明が無い)。未確認: 領域内の空洞の wake の横の蹴りを gs と比べていない、加速 + RAD | `C2b_sol_rf.sad`、C2 §2.2 |
| **スピン追跡 (POL + RAD、track6) の SOL 領域・UND・PHSROT** | v2 | ビームは 9 列 (7 列目 = atan(sz/sx)、8 列目 = sy)。スピンは RAD が ON のときだけ動く。通る要素: SOL 領域 (ずれ・回転・F1、領域内の DRIFT・QUAD・BEND・MULT・CAVI・SOL)、UND (磁場つきも。縦磁場と周期の途中で終わる端のポテンシャルの飛びの歳差も入る)、PHSROT (位相空間の写像だけでスピンは動かさない)。精度: 一様磁場の刻みは厳密、四極・多極の中は刻みの粗さ (1e-6〜5e-5。`EPS` を小さくすると細かくなる)。リングのスピンチューン: 平面リング (0.5 GeV) + 長さ 1 m のソレノイド (0・1・2.5 T) で cos(μ/2) = cos(π aγ) cos(φ/2) と 2e-5 | **Emittance[]・EMIT の POL の量 (SpinTune・Polarization など) は SOL 領域の中が K64 の式 (N6) のまま**: ずれて傾いたソレノイドではスピンチューンで 2.3e-4、n0 の向きで 4e-3 rad 違う。領域の外の Emittance[] のスピンチューンは N8 の分だけ追跡と違う (平面リングで 2.8e-5)。領域内の APERT、k64 モード | `C2c_spin.sad`、C2 §2.3 |
| **Polarization (POL) の平衡偏極度** (B55) | v2 | 平面リング 92.4% (Sokolov–Ternov)。誤差のあるリング: 脱偏極率が FLUC つきの追跡と 0.5〜1.6 倍 (3 例)。`emit6d.emittance6d(…, pol=True)['spin']` に `equpol`・`depol_rate`・`taudep`・`spin_spread`・`spin_orbit`・`equpol_k64` (K64 の式の参考値) | 線形 (SLIM 相当) の見積もり。高次の側帯共鳴、運動論的偏極 (b·γ∂n/∂γ の項) なし。エネルギー幅の大きいリングでは側帯の寄与が増えるはず。確かめたのは小さな 1 つのリング。POL のときは四極を細かく分割するので Emittance[] が遅くなる | A §3・§7 |
| **RADTAPER** | 両 (gs と同じ結果) | `RFSW; RAD; RADCOD; RADTAPER;` が ON のとき、Emittance[]・EMIT・TrackParticles・CALC6D が磁石の強さを閉軌道の局所の運動量に合わせる。`Twiss["DDP"]` に閉軌道の δ。閉軌道の z は 0、位相の基準は `EquilibriumPosition` (trf0)。gs の例: 軌道 −7.887e-5 → −2.198e-7 (pysad −2.173e-7: gs は収束判定が緩く、`CODCONV = 1e-4` で −2.17319735e-7 と 6 桁一致)、U0 8 桁一致、εx 1.7372e-7、εz 1.9077e-5 | CALC (4 次元の光学、K64 の qtwiss の radtaper は未翻訳)、SOL 領域内の要素、RFSW OFF。見つからないとき `???-RADTAPER-closed orbit with taper not found` を出してテーパーなしで続ける。追跡の閉軌道は Emittance[] と 1e-6 m 程度違う (K64 も同じ) | `C4_loss_and_taper.sad`、`C2d_taper_optics.sad`、C §2.5、C2 §2.4 |
| **LOSSMAP と損失** | 両 (gs と同じ結果) | `LOSSMAP;` で TrackParticles の戻りの flag の後に 2 行 (落ちた要素の l_track、落ちた周。生存は 0 と nend+1。POL のときは 10・11 行目)。POL + RAD の TrackParticles も LOSSAMPL・LOSSDZ で落ちる。gs の例と一致 (要素 {0,0,242,242,242}、周 {31,31,4,1,1}) | MAP 要素 (ExternalMap) が落とした粒子の要素番号 (0)、空間電荷が落とした粒子の要素番号。POL なしで高速追跡に乗らない要素の並びは SADError | `C4_loss_and_taper.sad`、C §2.4 |
| **MULT のスキュー成分** (`SK0 SK1 SK2 …`) | 両 (pysad の修正) | CALC・Emittance[]・追跡。gs と追跡 4 成分が 2e-13、Twiss (BX BY R1 R2 DX DY) が 1e-9 | — | §(d) N2、`C6_misc.sad` |
| **ResponseMatrix の TRPT** | 両 (gs と同じ結果) | 応答に sqrt(gammab(stt+1)/gammab(mon)) を掛ける。gs と全桁一致 (1 → 1.2 GeV) | — | `C6_misc.sad` |
| **RADLIGHT の軌道** (B65) | v2 | `RADLIGHTDIV = n;` で BEND の中の点数 (1/γ の角度あたり。既定 16)。3 GeV・0.26 rad の BEND で 2.5 万点 | SOL 領域を含む区間 (SADError)、COORD による幾何の変化、曲げのある MULT の中の点。DX・DY・F1・FRINGE のある BEND は中の点が 1e-6 m 程度ずれる | `C5_radlight.sad`、C §2.6 |
| **DELCOR (v2 の規格化)** | v2 | corlength > 0 の DX の rms = delx、DY の rms = dely (K64 は 1/√2)。CAVI・MARK にも入る | rms は多数回の平均としての値 (GCUT を掛ければ小さい)。相関長が格子間隔より短いと線形補間で rms が下がる。corlength=0 の枝は別の量 (分散 = delx²·s(1−s/周長))。`example/v2_C/C6_misc.pysad.txt` の DELCOR の行の値は √2 倍になる | `C2w_delcor_rms.sad`、C2_misc §2 |
| **旧式の WAKE コマンド** (B66) | v2 | `WAKE LONGITUDINAL 名前 s1 w1 … ENDWAKE` (横は `TRANSVERSE`)。ENDWAKE (または `;`) まで次の行に続く。CAVI だけ。名前はワイルドカード可。点数に上限なし (400 点で確認)。追跡は WakeFunction と同じ (TRPT と LWAKE/TWAKE) | 表の途中の `;`、CAVI 以外。ENDWAKE を書き忘れると次の `;` までを表として読んでエラー | `C2w_wake_table.sad`、C2_misc §1 |
| **RADINT の警告** (B74) | v2 | 場のある MULT・SOL があれば「入っていない」と表示 | — | A §1 |
| **EMITDIV・EMITDIVB・EMITDIVQ** | 両 (pysad の修正) | Emittance[] の放射の分割の倍率。MAIN の `EMITDIV = 20;` も効く。gs: `EMITDIV=20` で DampingRate が −.000286831320204 (pysad k64 と 12 桁一致) | — | §(d) |

### MEASURE・TRCOD と周辺

| 機能 | モード | 使える範囲と条件 | 使えない範囲 | 例・記録 |
|---|---|---|---|---|
| **MEASURE・TRCOD・LOGMEASURE・MAIN の `TRACK TURNS=-1`** (B41・B42) | v2 (K64 の gs は MEASURE が無出力、TRCOD が SIGSEGV) | 先頭が MARK のビームライン (輸送路、リング。リングは `CELL; CALC;` の後)。NP ≥ 2。MOMENTUM を設定する。輸送路は先頭 MARK に BX BY AX AY (と EX.. DX..) を書き EMITX EMITY DP SIGZ を MARK か FFS の変数で与える。**MARK に DP を書かないと DP は 0.01 (既定) で運動量の広がり 1% (一様) の分布**。追跡は TrackParticles と同じ経路 (RAD・FLUC・RFSW・APERT の損失)。`FFS["MEASURE 位置"]` は `{np, {重心 6}, {分散行列 6×6}}`。MEASURE の TRPT は追跡の間だけ ON (K64 は戻さない) | **gs との比較は無い**。`DRAW` の中の粒子追跡 (K64 でも死んだコード)。TRGAUSS OFF の分布と TWAKE/LWAKE ON の MEASURE は動くが独立な検証は薄い (LWAKE・TWAKE は一様・階段状のウェイクで重心と広がりを確認)。**強い結合のリング (モードが入れ替わる形) の MEASURE は未確認** (`Twiss["DETR"]` が常に 0: §(e))。MAIN の TRACK TURNS=-1 は粒子種の語 (PROTON など) を見ない | `C2m_*.sad`、C_measure、C2_measure |
| **TRCOD FIXCOD / MVCOD / STNCOD / ALIGN、ORBITJITTER、PUTMEASURE / GETMEASURE** | v2 | FIXCOD: EX EPX EY EPY DX DPX DY DPY を基準として保存し `Twiss["RDX",k]` (基準)・`Twiss["DDX",k]` (差) で読む (REFERENCE と同じ参照光学に置く。`REF;` は基準も上書きするので戻すには `TRCOD STNCOD;`)。ALIGN: 対象は BEND QUAD SEXT OCT DECA DODECA MULT の DX・DY だけ。ORBITJITTER: 今の Twiss の DX・DY に位置ごとの正規乱数。PUTMEASURE・GETMEASURE: 環境変数 SADPOUT・SADPIN のファイル (Windows は普通のファイル、Linux・WSL は FIFO を想定) | ALIGN の SOL・MONI。GETMEASURE の対象 (DX・DY) は判断。FIFO の待ち合わせは未確認。区間に BEND を含む ALIGN は平行移動と 4e-10 m まで | `C2m_trcod_reference.sad`、`C2m_putget_measure.sad`、C2_measure |

### 集団効果

| 機能 | モード | 使える範囲と条件 | 使えない範囲 | 例・記録 |
|---|---|---|---|---|
| **空間電荷 SPAC・PSPAC (+ SPCH 要素)** (B89〜B98) | v2 | `TrackParticles` の中。SPAC は `DRIFT` (`RADIUS>0`) と SOL 領域内のドリフト、MULT の `RADIUS` (ANGLE=0・VOLT=0・領域外。`SPACSTEP`)、PSPAC は `SPCH` 要素。電子 5 MeV〜50 MeV、1e8〜1e12 個/バンチ、ビーム半径 0.1〜2 mm、ダクト 3 mm〜3 cm。マクロ粒子は 1000 個以上 (チューンシフトを数 % で見るなら 4000 個以上)。条件: 粒子に z の広がり、`RFSW` ON (`NORFSW` はエラー。`RADIUS<0` の絞りだけなら走る)、ビームの外縁が 0.9 R 以内。FFS の変数 `PBUNCH`・`SPACSTEP`・`SPACNX/NY/NZ`・`PSPACNX/NY/DX/DY/NZ/DZ/NTURNCALC`・`WSPACBEAM`・`WSPACRADIUS`。フラグは FFS の最上位に書く (`Do[...]` の中では切り替わらない) | **k64 モードは空間電荷を追跡に入れない**。MULT の `ANGLE≠0`・`VOLT≠0`・SOL 領域内の MULT (メッセージを出して空間電荷なし)。四極・偏向・六極の中 (偏向のあるリングでは SPAC は全周に入らない。全周は PSPAC + SPCH か WSPAC)。光学・Emittance への空間電荷。MAIN の TRACK、DA、スピン追跡。複数バンチ、壁の抵抗、円でないダクト。未確認: 空間電荷とウェーク・放射の組み合わせ、バンチの静止系で相対論的になる場合 (γ·x' ≳ 1)、電子以外 (式は確かめたが SAD スクリプトからは使えない: `TrackParticles` が MASS・CHARGE を追跡に渡さない) | `C2s_mult_radius.sad`、`example/05_tracking/flag_SPAC` (作り直した)、C_spch、C2_spch |
| **WSPAC・SELFCOD** (強-弱のモデル) | v2 | 各要素の入口でガウス分布の力を蹴る。ビーム行列は既定で Emittance の平衡値、`WSPACBEAM = {εx, εy, σz[, σδ]}` で Twiss から作る。`WSPACRADIUS` で縦の力の電位の基準。電子 5 MeV (WSPACBEAM)・0.5 GeV (Emittance)、Δν = −3e-5〜−0.025。gs と 1% (0.5 GeV のリング、2e10 個) | 光学 (CALC)・Emittance への空間電荷 (ビーム行列は空間電荷なしの平衡値)。垂直のエミッタンスが 0 のリングは `WSPACBEAM` を与える。SOL 領域内・SELFCOD と閉軌道のあるリング・放射との組み合わせは値の確認なし | `C2s_wspac.sad`、C2_spch §2 |

### スクリプト・FFS・特殊関数 (モードによらない追加・修正)

| 機能 | モード | 使える範囲と条件 | 使えない範囲 | 例・記録 |
|---|---|---|---|---|
| **文脈** (`BeginPackage`・`Begin`・`End`・`EndPackage`) | 両 (k64 も名前空間を作る) | K64 (init.n と gs の実測) の規則: 名前は「現在の文脈 → `$ContextPath`」の順に探し、無ければ現在の文脈に作る。`$Context`・`$ContextPath`・`Context[sym]`。大域に公開する関数は `` `Name[…] := … ``。組み込み関数・FFS の関数はそのまま使える。`Names["Stat`*"]` は v2 が文脈の記号の一覧、k64 が gs と同じ `{"Stat`"}`。SAD 自身の Packages は従来どおり文脈を無視して大域で読む。K64 との差 (意図的): pysad の System` は K64 より広い、表示は文脈名を `$ContextPath` にあるものだけ省く、`Symbol["Stat`x"]` は記号を返す | `$Context = …`・`$ContextPath = …` への代入、`WithContext`。入れ子の文脈は確かめていない | `B2_package_context.sad`、B §2.4 |
| **Fork・Shared・Wait・Kill・Sleep** | 両 (Windows の Fork は v2 でエラー、k64 は従来どおり pid を返す) | `OpenShared`・`Shared[s]` (数・文字列・記号・リスト・式・規則)・`SharedSize` (K64 と同じ式。`Shared::toolarge`)。POSIX (Linux・WSL・macOS) の `Fork[]` は `os.fork` で分ける。`Wait[]`・`Wait4`・`Kill[pid]` (この pysad が作った子だけ)・`Sleep[t]`。WSL で example と 3 個の子の試験 | Windows の Fork (`System::error` で打ち切り。並列は `Parallelize`、pysad は逐次)。ラティスを使う重い計算を子で走らせる使い方 (numba のスレッドと fork の併用) は未確認。子の `Print` は出力先が StringIO のとき親に届かない | `B3_shared_system.sad`、B §2.5 |
| **System** | v2 (許可なし) | 既定では実行しない (`$Failed` を返し案内を出す)。許可: 環境変数 `PYSAD_ALLOW_SYSTEM=1` (起動前)、または `pysad.miscfn.allow_system(True)`。SADScript からは変えられない。許可があるとシェルで実行し終了コード (成功 0) を返す。`[ -e path ]`・`[ -d path ]`・`[ -f path ]`・`[ -r path ]`・`[ -s path ]`・`test …` は許可なしでも判定する。`PYSAD_SYSTEM_SHELL` で sh・bash を指定 | k64 モードは許可なしで 0 (「実行しません」)、許可ありで戻り値 1。SAD の Packages が呼ぶ個々のコマンド (`import`、`latex` など) は未確認 | `B3_shared_system.sad`、B §2.6 |
| **複素引数** (Bessel J/Y/I/K (実数の次数)、Hypergeometric0F1・1F1・2F1 (実数のパラメータ)、Airy) | 両 | 積分表示・級数・Wronskian と 1e-12〜1e-13 | 複素の PolyLog・HypergeometricU・Zeta (scipy に実装が無い・分岐の向きを確かめる参照が無い)。`PolyLog` は非整数次で負の s が大きい (s < −10 程度) 精度は未確認 | `B1_special_functions.sad`、B §2.2 |
| **複素の FindRoot** (B61) | v2 | 式が変数について正則と仮定した Newton 法 (前進差分のヤコビアン + 最小二乗、減速つき)。連立可。`MaxIterations`・`AccuracyGoal` が効く | 変数の範囲 `{var, ini, {min, max}}` は複素では無視。Abs・Conjugate・Re・Im を含む式は収束を保証しない。実数の初期値で複素の値になる式は従来どおり `General::wrongval` | B §2.3 |
| **MOMENTUM の警告** (B67)、**`MON`・`MCOD`・`TCOD` の警告** (B22) | v2 | MOMENTUM 未設定で運動量を使う計算に入ると 1 度だけ警告 (CALC と幾何では出さない。SADScript の式の中の `MOMENTUM` の直接の使用では出ない) | — | `B4_ffs_v2.sad` |
| **MAIN の代入と `FFS` の文の読み方** | 両 (pysad の修正) | gs の動作 (小さな台本 15 本) に合わせた: MAIN の `NAME = 式;` (数 [単位]、`+ - * /`、`^` (右結合)、関数 `SIN COS TAN ATAN SQRT EXP LOG LN COSH`、`PI`、`MOMENTUM`) が要素定義と FFS に届く。`FFS USE RING;`、行の途中の `FFS`、`FFS;NPARA=3;` の読み。32 個の式を gs の値と 2e-14 | gs で構文エラーになる形は代入しない (`Sqrt[16]`、`ABS(…)`、`2*-A` など。`(1+2)*3` だけは読む)。残した違い: `Ffs;` (大小まじり)・`ffs use ring;` は pysad は読む、`FFS` の次の行が USE 以外のとき、gs が FFS の終了後に MAIN の文として読む残りは pysad は読まない、`A = 3 B = 4;` は代入しない、`LINE[{"EMITY",…},1]` は MARK の値 (0) を返す、要素定義の中の MAIN の関数 (`L=2*SIN(A)`) は未確認 | `C2w_main_script.sad`、C2_misc §3・§4 |
| **NOKEEPEXP・RENUMBER の INDEX・`LINE["INDEX",…]`** | 両 (gs と同じ結果) | `NOKEEPEXP` のとき EXPAND/NOEXPAND を付けない CALC・GO は `calexp = .not. tfvcomp()`。INDEX 属性が 2 以上の要素の番号はその値から。gs の実測と一致 | `tfvcomp` は常に偽 (成分ごとの FREE が無い)。RENUMBER 後の `LINE["POSITION",…]` の gs の癖は写していない | `B4_ffs_v2.sad`、B §1 |

## (d) pysad v1 の不具合の修正

K64 のバグではなく、pysad 自身の誤り・未対応だったもの (K64 の番号 B72 以降には入れていない: [`21_sad_bug_report.md`](21_sad_bug_report.md) §6)。記載は「モードによらない」ものを含む。

| 項目 | v1 の不具合 | 修正 | 根拠・試験 |
|---|---|---|---|
| `LegendreP[n,m,z]` の奇数 m の符号 | scipy の `lpmv` に (−1)^m を掛けていた (`lpmv` は Condon–Shortley の位相つき。`LegendreP[3,1,.4]` が −.2749…) | 掛けない。**+.2749545416973507** | 式の m = 1 ± 1e-6 の値と閉じた式 (`test_b45_legendre_integer_m_is_limit_of_formula`)。gs は NaN (B45) で参照が無く気づけなかった。モードによらない |
| MULT の `SK0 SK1 SK2 …` (N2) | `lattice._mult_params` が `KsL` を読まず、スキュー成分を黙って無視 (CALC・Emittance[]・追跡) | `akr(n) = cmplx(K_n, SK_n)·cr1^(n+1)`、cr1 は複素 K1 の主軸 (tpara) | gs: 追跡 4 成分が 2e-13、Twiss が 1e-9 (`test_v2_C_sol.py::test_mult_skew_components_match_gs`)。モードによらない |
| `CalculateOptics` の TRPT の加速 (N3) | CALC は入れるが CalculateOptics は入れない | TRPT かつ RFSW の加速を区間に入れる (区間の途中からでも) | gs: 終端の BX AX BY AY EX EPX NX NY が全桁、CalculateOptics − CALC = 0 (`test_v2_C2_optics.py::test_calculate_optics_includes_trpt_acceleration`) |
| 加速区間の分散 (N9) | qmat2twiss の gr = gammab(出口)/gammab(入口) の向きが逆 (入口に分散のある加速区間の EX EPX が gs と合わなかった: 0.1109 / 0.0734) | K64 のとおり gammab(出口)/gammab(入口) を入口の分散に掛ける | gs と全桁 (薄い空洞・厚い空洞)。`test_dispersion_through_accelerating_cavity_matches_gs` |
| EMITDIV・EMITDIVB・EMITDIVQ | pysad は無視していた | Emittance[] に効く (MAIN の `EMITDIV = 20;` も) | `EMITDIV=20` で DampingRate が gs と 12 桁 (`test_emitdiv_is_used`、`…_in_main_section`) |
| SOL の L と S 座標 (B63) | pysad は S 座標 (LENG) にだけ SOL の L を入れていた (46.9) | 入れない (gs の LENG 46.6)。SOL に L のキーワードは無い (SADHelp) | `test_v2_A_ffs.py::test_sol_has_no_length`。モードによらない |
| 特殊関数 | 非整数次の `PolyLog` が z < −12・z > 60 で精度が落ちる (z = −140、s = −1.5 で相対誤差 1。\|z\| ≥ 150 は NaN)。複素引数の Bessel・超幾何・Airy が TypeError | `_polylog_far` (反転公式・積分表示)。複素引数の対応 (§(c)) | s ≤ 7.9、\|z\| ≤ 1e12 で相対 1e-11 |
| `SharedSize`・`LINE["INDEX",…]`・INDEX 属性の番号・`h::tag` の `General::tag` への継承・POSIX の `Fork`・`Sleep` | `SharedSize` が 1 (gs は 32)、`LINE["INDEX",…]` が 0、INDEX 属性が効かない、`Fork` が分けない、`Sleep` 未定義 | gs の動作に合わせた (§(c)) | gs の実測 |
| `fn_track` (B33 の修正) | POL + RAD のとき後ろの `_compat.v2()` が UnboundLocalError (関数内の import が原因。k64 の `test_spin_track.py` が 5 件落ちた) | モジュール先頭の import だけにした (`measure.py`・`_exec_word` も同様) | 試験 |
| `Element["RADIUS", ..]` (DRIFT・MULT) | 0 を返していた | 定義の値を返す (gs と同じ。追跡・光学の値は変わらない) | `test_mult_radius_readable`、`test_element_radius_readable`。モードによらない |
| 描画の試験 (`test_plot.py`) | Tk がマウスの下の項目に付ける `current` タグで判定が変わる | 読み取りで `current` を外す (判定は同じ) | — |
| 例 `example/05_tracking/flag_SPAC` | `NOSPAC;`・`SPAC;` を `Do[...]` の中に書き、フラグが切り替わらず「gs でも SPAC の効果が出ない」と誤記 | 最上位に書き、解析式と比べる例に作り直した (gs は力が 1/γ: B89) | — |

## (e) 要判断・未解決

各記録の「要判断」「残した違い」「使えない範囲」を、解決済みのものを除いて集めた。「決定済み」は親の判断で変更していないもの。

### 放射・エミッタンス・スピン (記録 A)

- **スピン追跡の 500 周後の `spin angle` (`example/07_spin_solenoid_special/flag_POL_RADPOL`)**: v2 の既定では 1.1516 rad で、gs と k64 モードの 0.8794 rad と違う。原因は特定していない (A の B58・B56・K1 つき BEND の追跡の修正のいずれかによる意図した違いの可能性が高いが、独立な根拠で 500 周後の値を確かめていない)。小さなリングの 1 周あたりのスピンの回転角を解析式 (スピンチューン aγ と、ソレノイド・K1 つき BEND の端の項) で確かめる試験が要る。
u_s=a\gamma$ とソレノイド・K1 つき BEND の端の項) で確かめる試験が要る。

| 項目 | 見つけたこと | 直さなかった理由・状態 |
|---|---|---|
| tradke の dbeam の b = bzh·0.5 (B106) | 力学的運動量は px + bzh·y なので b = bzh のはず | bzh の定義の別の約束がある可能性を排除できない。SOL 内の放射のビーム行列の「ビームの広がり」の項だけに効く |
| tbgrote の trans(5,6) (B109) | コメントアウトされた微分の式ではなく別の式。dr=0 でも dpzi·(pr+1)/pzi が残る | 作者が意図して置き換えた式。2 次の小さい量 |
| RADINT の CKL1Y (B107) | sin の誤記に見える (結合のあるときの I4V) | 結合のあるリングの独立な I4V の計算を持っていない |
| RADINT の BEND の K0 (B74) | 無視される (ANGLE だけ) | 仕様かもしれない |
| RADINT の空洞の位相 | VOLT を単純に足す (PHI を見ない) | 位相の違う空洞のベクトル和にするかは仕様の判断 |
| K64 の spdepol (2 次の摂動、シンクロトロン側帯) | 大きさが追跡と 3〜5 桁合わない (側帯の項に 2π が余分、1/減衰率の因子)。v2 は線形の式に置き換えた | K64 の式を直すには導出が要る |
| v2 の Polarization の限界 | 軌道の振幅について 1 次 (SLIM 相当)。追跡との一致は 3 例で 0.5〜1.6 倍 (16 粒子、統計 ±25%)。高次の側帯共鳴、運動論的偏極なし | 確かめたのは小さな 1 つのリングだけ |
| CALC (tbende) と追跡 (tbend) の、ずれた BEND | 矩形磁石を DX=2 mm ずらすと CALC のチューンが 2.9e-5 変わる (追跡は変わらない。追跡が正しい) | tbende の 2 次の幾何の項。場所を特定していない |
| K1=0 の BEND の追跡のスピン (tbend + tradk) | 1 個の BEND のスピン回転角が Emittance[]・tbendi (v2) と 2.2e-5 rad 違う | 放射の分割が粗い K64 の精度の範囲 |
| B21 の K64 のソース | tcod.f の rfsw=.false. の枝は δ を動かさない。z を閉じる δ は temit の呼び出し側で決まるはず | 呼び出し側を追っていない。pysad の動作は gs の実測と合い、物理としても妥当 |
| ANGLE=0 の BEND・K0 つき MULT の端の角 (B108)、tthine の bsir0 (B110) | B75 と同じ作りのはず / 次元が合わない | 独立な値が無い / 影響なし (直していない: ANGLE≠0 の BEND だけ B75 を直した) |
| 小さな Jx の残り | 既定の分割 (EMITDIV=1) では Emittance[] の Jx が 1.5e-3 大きい (K64 の精度)。Jx+Jy+Je = 3.9990 (U0/2E の 2 次) | `EMITDIV=20;` で 1e-4 まで合う |
| 範囲の限界 | 大きなリング・挿入光源・強い結合では確かめていない。RADINT の V 面 (I4V) は確かめていない。SOL 領域内の MULT・BEND のスピンは (Emittance[]) 未対応 | — |

### 追跡・要素 (記録 C・C2)

| 項目 | 状態 |
|---|---|
| B17 (NaN の返し方) | 決定済み: 現状のまま。asinz・gfloor・sqrt1 の下限は、発散した粒子が落ちるまでの座標の埋め方にだけ効く。「NaN になった粒子をその場で落とす」に変えると落ちる要素・周が K64 と変わる (LOSSMAP の値)。マッチングの fuzz (tffsmatch) は未確認 |
| DELCOR の CAVI の符号 (B40) | 決定済み: v2 の +delx (K64 は設計の DX − delx)。SOL の式を写した意図があるなら逆 |
| DELCOR の corlength≠0 の規格化 (B99) | v2 は rms = delx に直した (作者の文書による裏づけは無い。推定)。K64 の値が必要なら k64 モードか delx に 1/√2 を掛ける |
| UND の光学 | 決定済み: v2 は磁場なしでも光学にドリフトとして入れる (追跡と揃える)。K64 の「光学では恒等」を仕様とみるなら磁場つきのときだけにする選択もある |
| B33 の落ちた要素 | v2 は「始点の 1 つ前の要素 (周回で数える)」。K64 の 2 回目以降は「前の追跡の最後の要素」 |
| スピン追跡の損失の k64 の最初の呼び出し | K64 の B33 の癖 (落とす印 0 の粒子が後で入れ替えで戻る) は track6 の経路では近似 (gs 要素 60、pysad 15)。高速追跡 (track6n) は全桁一致 |
| 未対応のまま (SOL・UND まわり) | スピン追跡の経路 (track6) での TRPT の加速 (RF は入るが基準運動量は変えない)、領域内の APERT のスピン追跡、UND の端の磁極の整形、CALC (4 次元) の RADTAPER (K64 も無し)、bz = 0 の領域と領域の外の RF つき MULT の Emittance[] (C2 §2.2。C2 §1 の項目 2d は「bz = 0 の領域の中の CAVI・MULT は領域の外と同じ写像」とも書いており、**記録内の表現に食い違いがある: 要確認**) |
| Emittance[]・EMIT の POL の量 (SOL 領域内) | K64 の式 (N6 の bsi) のまま。`spin.tradke_spin` の `xpzb = xpz + (bsir0 + bzh·2·al)·pr²` を直すには 6 次元の微分 dxpzb も一緒に直す必要がある (未着手) |
| N12: `Twiss["DETR"]` が CALC の後で常に 0 | **pysad の未修正の不具合** (`ffscmd._rows_from_tws` が detr を渡していない。k64・v2 とも)。結合の強いリングの MEASURE に効く。共有ファイルの CALC の中心なので直していない |
| N13 (trots の入口の DZ)、AUTOPHI の意味 (N11) | K64 のソースの読み。gs で確かめられたのは TRPT の AUTOPHI だけ。リング (TRPT でない) の追跡の AUTOPHI は K64 と違う |
| RADTAPER の追跡の閉軌道 | Emittance[] の閉軌道から追跡を始めても、1 周で x が 4e-7 m、z が 2e-6 m 動く (テーパーなし)。K64 の追跡と tturne の放射の入れ方の違い (pysad の k64 の追跡は gs と 1e-12 で合う)。直していない |
| 領域内の未確認 | 領域内の MULT・CAVI の POL、RF つき MULT の RAD、bz ≥ 3 1/m を超える場、領域内の空洞の wake の横の蹴り (gs と未比較)、加速 + RAD。gs は SOL + 四極 + RAD で保存値と食い違う (B68) |
| gs との差の未解決 (docs/46) | BEND の `F1 FRINGE=3` + RAD の追跡、L≠0 の CAVI + `DX`・`DY`、EMIT のあとに `USE` したときの追跡 は v2 でも未解決 |

### MEASURE・TRCOD (記録 C_measure・C2_measure)

| 項目 | 状態 |
|---|---|
| gs との比較 | 1 件も無い (gs は MEASURE が無出力、TRCOD が SIGSEGV)。実装は K64 のソースの読み。乱数を引く順はソースどおりだが粒子ごとの値を gs と比べたことはない |
| MEASURE の TRPT の副作用 | 決定済み: v2 は追跡の間だけ ON にして戻す (K64 は戻さない)。K64 のとおりにしたければ MEASURE の後に `TRPT;` と書く (または `cmd_measure` の `finally` の 1 行を外す) |
| QUIET の重心合わせ (B87)・TWAKE/LWAKE の z = 0 (B88) | v2 は直した (K64 の式の誤りと判断。名残と判断)。TRCOD の途中の統計 (K64 は装置番号 0) は出さない |
| GETMEASURE の対象 (B103) | DX・DY の条件と判断。旧い DDX.. (基準との差) だった可能性があり、その場合の目標値は基準の分だけ違う |
| FIXCOD の基準の置き場 (B100) | REFERENCE と同じ参照光学に置いた (`measure._ref_rows` の 1 か所で替えられる)。ALIGN の対象に DECA DODECA MULT を足した (`measure._ALIGN_V2` を `_ALIGN_K64` にすれば K64 の BEND QUAD SEXT OCT だけ) |
| 前任の試験の古くなった 1 行 | C2 の記録では `test_v2_C_measure.py::test_measure_particle_loss_and_errors` の `'?TRCOD FIXCOD: not supported'` の 1 行 (FIXCOD の実装で出なくなる) と `test_v2_C_delcor.py` の ⟨\|c\|²⟩ の期待値が古くなっていた。現在の試験は `'  trcodf'` と 2.0 を見ていて、全試験は通る (解決済み) |
| 未翻訳 | `DRAW` の中の粒子追跡 (K64 でも死んだコード)、PUTMEASURE・GETMEASURE を FIFO で使う待ち合わせ (未確認)、ORBITJITTER のパイプ |

### 集団効果 (記録 C・C_spch・C2_spch)

| 項目 | 状態 |
|---|---|
| 空間電荷の k64 モード | 決定済み: K64 の動作 = 空間電荷なし。PSPAC だけ v2 の実装を k64 でも使うかは未決定 (gs に 3〜4 桁で近づくが「k64 = 空間電荷なし」が変わる) |
| `TrackParticles` が `MASS`・`CHARGE` を追跡に渡さない | 電子以外は追えない (空間電荷と関係のない pysad の既存の制限。他の追跡にも影響するので親の判断) |
| 光学・Emittance への WSPAC (K64 の `twspace` の反復) | 未対応。入れるとビーム行列が自己無撞着になり Twiss にも空間電荷のチューンシフトが出る (`emit6d`・`calc6d` の変更が要る) |
| 空間電荷の未確認 | ウェーク・放射との組み合わせ、SOL 領域内の四極・MULT、`ANGLE≠0` の MULT、静止系で相対論的な運動、放射 + 空間電荷の MULT (走るが PBUNCH=0 の一致だけ)。PSPAC の gs との 16% の差 (直線の例で x と y が違う) は原因を調べていない |
| 文書の更新 | `docs/46`・`docs/01` などの「SPAC は読むだけ」「gs でも効果が確認できない」の記述 (別の担当が更新) |

### 特殊関数・スクリプト・FFS (記録 B・C2_misc)

| 項目 | 判断 | 理由・調べたこと |
|---|---|---|
| B24 `LINE["LENGTH",loc]` | 値は返さず案内だけ足した | K64 のソースが `LENGTH` の引数を 1 個に限る。位置を付けたときの意味 (位置か、軌道長か) がソースからもヘルプからも決まらない |
| B22 MON・MCOD・TCOD の機能 | 実装せず警告して読み飛ばす | K64 にコマンドが無く、古い SAD の MON の書式・意味を確かめられる資料が無い。軌道補正は CorrectOrbit パッケージ (docs/34) |
| `Out[n]:=` の gs の表示 | 写さない | gs は未定義の記号 1 つの文で前の `Out[n]` を繰り返す。K64 の表示の癖 |
| NOKEEPEXP の `tfvcomp` | 常に偽 | pysad の `FREE` はファミリー単位で成分ごとの変数が無い。成分ごとの FREE を実装するときに `FFS._tfvcomp` を直す |
| RENUMBER 後の名前の引き方 | pysad は `LINE["NAME"]` の表示と一致する位置を返す | gs は RENUMBER のあと `LINE["POSITION","QF.1"]` が `LINE["NAME"]` の表示と食い違う (内部の表の更新の癖。INDEX の有無によらず以前からの差) |
| 複素引数の PolyLog・HypergeometricU・Zeta | 未対応 (SADError・TypeError) | scipy に複素の実装が無い、分岐の向きを確かめる参照が無い |
| `Parallelize` の本当の並列 | 逐次のまま | POSIX では Fork で実装できるが、numba のスレッドと fork の併用を確かめていない |
| `Symbol["x"]`・`Message[f::usage]` の gs との差 | そのまま | gs は `Symbol` が値を返し、`Message[f::usage]` が `???f::usage:  本文` を表示する。pysad は記号を返し、何も表示しない |
| 文脈の未確認 | `$Context = …` への代入、`WithContext`、入れ子の文脈 | §(c) |
| `Ffs;`・`ffs use ring;` など | pysad が gs より寛容 | §(c) の「残した違い」 |

### K64 側の未解決 (docs/21)

B68 (gs の 10-19 ビルドが SOL + 四極 + RAD で保存した参照と食い違う: 現在の gs が領域内の四極の力を入れていない疑い)、B60 (CSROYMatrix が gs で wrongtype)、B49・B51 (RadiationField・RenumberElement が gs で SIGSEGV) などは v2 でも比べられない。B72〜B110 のうち確度が「疑い」のもの (B88、B99、B106〜B109) は K64 の作者の意図を確かめていない。

## (f) 試験と回帰の状態 (2026-11-17 20:30 に実行)

| 項目 | 結果 |
|---|---|
| 全試験 | `py -3.12 -m pytest pysad/tests -q` (旧式の 3 本 `test_bend_feeddown.py`・`test_coord.py`・`test_sext_radiation.py` を `--ignore`、環境変数 `PYSAD_SAD_PACKAGES` つき、互換モードは k64 が既定): **937 通過・4 skip・失敗 0** (941 件収集) |
| v2 の試験 | `test_v2_A_*.py` 67 件、`test_v2_B_*.py` 52 件、`test_v2_C_*.py` 114 件、`test_v2_C2_*.py` 130 件 (計 363 件、29 ファイル)。v2 の動作は `compat.set_mode('v2')` で明示して試験する |
| 回帰 | `run_regression.py` (`QUICK=1`): **50 ベンチ FAIL 0**。k64 モードで走る (gs との比較) |
| skip | GUI 系は `PYSAD_SAD_PACKAGES` が無いと skip (環境変数を付けると走る)。`PYSAD_SLOW=1` で遅い試験が走る。4 skip の内訳は未確認 |
| 条件 | Windows (Python 3.12)。`Fork` の試験は WSL (Ubuntu、Python 3.14、numba は代用の空モジュール) で確認 |
| 限界 | 試験の参照は解析式・独立な数値計算・gs。v2 の値は gs では確かめていない (gs が K64 のバグを持つため)。MEASURE・TRCOD は gs との比較が無い。実機の大きなリング・挿入光源・強い結合では未確認 |

## 記録と例

| 記録 (`history/`) | 内容 |
|---|---|
| [`v2_changes_A.md`](history/v2_changes_A.md) | 放射・エミッタンス・スピン・RADINT・tchge (B06〜B16、B39、B55・B56・B58・B64 ほか、Jx の 2.4% の差の結論、平衡偏極度) |
| [`v2_changes_B.md`](history/v2_changes_B.md) | 特殊関数・数値、SADScript、FFS の小項目 (文脈、Fork・Shared、System、複素の FindRoot) |
| [`v2_changes_C.md`](history/v2_changes_C.md)、[`C_measure`](history/v2_changes_C_measure.md)、[`C_spch`](history/v2_changes_C_spch.md) | 追跡・集団効果・要素 (UND、PHSROT、SOL 領域内、LOSSMAP、RADTAPER、RADLIGHT、旧 WAKE)、MEASURE・TRCOD、空間電荷 |
| [`v2_changes_C2.md`](history/v2_changes_C2.md)、[`C2_measure`](history/v2_changes_C2_measure.md)、[`C2_misc`](history/v2_changes_C2_misc.md)、[`C2_spch`](history/v2_changes_C2_spch.md) | 残りの仕上げ (ずれたソレノイド、領域内の空洞、スピン追跡、RADTAPER の収束、TRCOD の残り、旧 WAKE の表・DELCOR・台本の読み方、MULT の RADIUS・WSPAC) |

| 例 | 内容 |
|---|---|
| [`example/v2_B/`](../example/v2_B/README.md) | `B1_special_functions`・`B2_package_context`・`B3_shared_system`・`B4_ffs_v2` (特殊関数、文脈、共有メモリ・System・Fork、FFS の小項目) |
| [`example/v2_C/`](../example/v2_C/README.md) | `C1_undulator`、`C2_phsrot`、`C3_sol_region`、`C4_loss_and_taper`、`C5_radlight`、`C6_misc`、`C2a_sol_tilt`、`C2b_sol_rf`、`C2c_spin`、`C2d_taper_optics`、`C2s_mult_radius`、`C2s_wspac`、`C2m_trcod_reference`、`C2m_putget_measure`、`C2w_wake_table`、`C2w_delcor_rms`、`C2w_main_script` |
