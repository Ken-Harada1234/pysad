# v2.0.0 の変更記録: 担当 C2 (担当 C が残した項目の仕上げ)

作成: 2026-11-17。担当 C の記録 ([v2_changes_C.md](v2_changes_C.md)) の §6「残っているもの」と §3「要判断」を仕上げた分です。
分けて書いたもの: 空間電荷 (MULT の RADIUS・WSPAC・例 flag_SPAC) は [v2_changes_C2_spch.md](v2_changes_C2_spch.md)、
MEASURE 系 (TRCOD FIXCOD/MVCOD/STNCOD/ALIGN、PUTMEASURE・GETMEASURE・ORBITJITTER) は [v2_changes_C2_measure.md](v2_changes_C2_measure.md)、
旧式の WAKE の表の行またぎ・DELCOR の規格化・台本の読み方 (`FFS;NPARA=1;`、MAIN の代入) は [v2_changes_C2_misc.md](v2_changes_C2_misc.md)。

- 互換モード: `pysad.compat` (`PYSAD_COMPAT=v2|k64`)。「k64 モード」の列が「同じ」のものは切り替えていない (pysad の不具合の修正・追加)。
- 根拠: 解析式、独立な数値積分 (scipy の solve_ivp: ローレンツ力、Thomas-BMT、縦電場の中の加速)、座標変換の合成、極限の一致。
  gs (SAD 1.4.4.2k64、2026-10-19 ビルド、WSL) は、K64 が正しく動く部分 (TRPT の光学、RADTAPER、台本の読み方) の参照と、K64 の動作の確認に使った。
- 試験: `pysad/tests/test_v2_C2_*.py` (件数は §5)。走る例: `example/v2_C/C2*.sad` と出力 (v2 モード)。
- 使った条件: 電子 (追跡の一部の試験は電荷 +1 の既定値)。0.5・1・3 GeV、周長 91 m の FODO リングと 4〜6 m の輸送路、
  ソレノイド 1〜3 T (bz = 0.3〜3 1/m)、空洞 20〜200 MV (1.3 GHz)・4 MV (500 MHz)。実機の大きなリングでは確かめていない。

## 1. 一覧

| 項目 | K64 の動作 | v2 の動作 | 根拠 (式・数値) | 試験 | k64 モード |
|---|---|---|---|---|---|
| 1a. SOL のずれ・回転 (`DX DY DZ CHI1..3`、tsol.f の trots) | tsgeo が、GEO≠0 の SOL の `DX DY` (位置)・`DPX DPY` (向き) の軌道を領域の中で追い、両端の SOL の CHI1..3 と、もう一方の SOL の DX DY DZ を決める。追跡は入口と出口で trots (回転 + 新しい基準面までの直進 + ずれ)。pysad v1 は幾何だけ対応で、追跡・Twiss・Emittance[] は NotImplementedError | 追跡 (高速追跡と track6)・CALC の Twiss と軌道・Emittance[]・EMIT に入る。入口の z にも −dv·DZ を入れる (K64 は出口だけ) | 実験室系の合成: 幾何 (GX GY GZ と軸) からソレノイドの面と軸を置き、[直進 → 端の蹴り → 一様磁場 (+ 四極・双極) の中のローレンツ力の積分 → 蹴り → 直進]。ドリフトだけの領域で 2e-12 (x px y py z)、四極・BEND・MULT つきは 1e-8 (カーネルの運動学的な高次の項)。GEO が入口側・出口側の両方、30 MeV (v≠c の z の項)。設計軌道の粒子は出口で 0 (1e-13)。BZ = 0 なら入口と出口の変換が打ち消して直進と同じ (1e-14) | `test_v2_C2_sol.py::test_tilted_offset_solenoid_matches_lab_frame_lorentz_integration` (5)、`…field_free_limit`、`…with_bend_and_mult_inside` | 従来どおり (NotImplementedError) |
| 1b. SOL の `F1` (端の放射、tradks) | 端の前後の力学的運動量の向きの変化を、長さ F1 の一様な偏向として tradk に渡す (RAD かつ F1≠0 かつ bz が変わる端)。pysad v1 は NotImplementedError | 同じ式で実装 (追跡・Emittance[] の tradke) | Δδ = −(2/3) r_e γ³ (Δbz·r/2)²/F1。入口 (1e-5)、領域内の SOL・出口 (ラーモア系で求めた半径で 1e-4)、F1 を 2 倍で半分。リングの U0 の増分 (2.11 eV) が式と 2e-3 | `test_sol_f1_edge_radiation_matches_uniform_kick_formula`、`test_ring_with_tilted_solenoid_calc_orbit_tunes_and_emittance` | 従来どおり |
| 1c. 領域の幾何 (tsgeo の BEND・MULT・CAVI、K1=0 の四極)、幾何の fit (tdgeo) | tsgeo1 は領域内の BEND・MULT・CAVI も軌道を追う。tdgeo (fit の微分) は領域内の変数にも領域の外と同じ式 (軸方向の伸びだけ) を使う。pysad v1 は領域内の BEND・MULT・CAVI と K1=0 の四極で幾何が作れなかった | 領域内の BEND・MULT・CAVI・K1=0 の四極も軌道を追う (追跡と同じカーネル)。fit の幾何の条件は、変数が領域の中にあると中心差分 (軌道の回転の変化が入る) | 設計軌道の粒子が出口で 0 に戻る (領域内に BEND + MULT、1e-12)。ほかの粒子はローレンツ力の積分と 1e-8。fit (GY・GZ) が領域内の四極の K1・ドリフトの L で収束し、目標に 1e-8 | `test_tilted_solenoid_with_bend_and_mult_inside`、`test_geometry_fit_through_tilted_solenoid_region` | 従来どおり |
| 1d. リング (CALC・Emittance[]) | (K64 は対応) | ずれて傾いたソレノイドのあるリングで CALC・Emittance[] が走る | CALC の閉軌道は領域の外で 0、中はソレノイドの系で (DX, DY) から。CALC のチューン = 追跡の 1 周の行列の固有値 (2e-7)。減衰率の和 = −2 U0/E (2e-3)。BZ → 0 なら SOL なしのリングと同じ (チューン 1e-6、ε 1e-5) | `test_ring_with_tilted_solenoid_…`、`test_ring_tilted_solenoid_zero_field_limit` | 従来どおり |
| 2a. 領域内の CAVI・RF つき MULT の TRPT の加速 | tmultiacc を bz つきで呼ぶ。bz は要素ごとに、その入口の基準運動量で規格化 (tfbzs)。pysad (担当 C) は明確なエラー | 実装 (追跡: tmultiacc と同じ重みの分割・出口で基準運動量を付け替え。CALC: 同じカーネルの数値微分)。加速の後の領域の要素は bz·p0(入口)/p0(その場所) | 一様な縦電場 Ez = −(V/L) sin(ω(t − t_ref) + φ) と Bz の中の運動方程式を基準粒子と一緒に積分: 出口の運動量 1e-14、座標 5e-10 (20 MV)、3e-8 (200 MV: K64 の空洞の分割の精度)。CALC (TRPT) の軌道 = 追跡 (2e-9) | `test_v2_C2_solrf.py::test_acceleration_in_solenoid_matches_equation_of_motion` (2)、`test_trpt_calc_orbit_through_accelerating_region_equals_tracking` | 従来どおり (領域内の CAVI は未対応) |
| 2b. 領域内の CAVI・RF つき MULT の wake | tsol.f の ewak: CAVI と VOLT≠0 の MULT は tmultiacc の中で分割ごとに txwake | 実装 | 縦の wake によるエネルギーの変化の合計が、領域の外の同じ空洞と同じ (bz → 0 で 1e-6、BZ = 1 T でも同じ)。横の蹴りの合計は 2% 以内 (分割して入るので位置の変化は少し違う) | `test_wake_inside_region_cavity` (2) | 従来どおり |
| N10 (新規、K64) wake の重み | tcav.f・tmulti.f: `fw=fw0*(wsn+ws(n-1)*.5d0)`。空洞が 2 つ以上に分割されると wake の蹴りの合計が (1.5n − 0.5)/n 倍 | 台形の重み fw0·(ws(n) + ws(n−1))/2 (合計 1)。領域の外の空洞・RF つき MULT も | 位置によらない表では、エネルギーの変化が空洞の電圧 (= 分割数) によらない: 20 MV (1 分割) と 200 MV (4 分割以上) で 5e-3。k64 は 1.25〜1.5 倍 | `test_wake_weight_sum_is_one_in_v2_and_one_and_half_in_k64` | 再現 |
| 2c. AUTOPHI | 追跡 (tturn.f・tcav.f): `ph = DPHI + DZ·w` を足し、`phis = merge(0, merge(phic, w·trf0, autophi), trpt)`。TRPT では空洞の位相の基準が設計軌道の到着時刻 (Twiss の DZ) になる。TRPT でないと phis = phic で位相の差が消える (PHI・DPHI・DZ が効かない)。光学・Emittance (tcave.f・tmulte.f): 軌道は位相 phic で通る (phii = phic)。SADHelp に説明なし。pysad は領域の内外とも未対応 (NotImplementedError) | 「位相の基準 = 設計軌道の到着時刻」として、追跡 (CAVI・RF つき MULT、領域の内外、TRPT でもリングでも phis = 0)・CALC (TRPT の加速)・Emittance[] に実装 | 設計軌道が DZ = 3 mm 早く着く輸送路: AUTOPHI なしは DDP = 5.2e-4 (CALC = 追跡)、ありは設計軌道の粒子の δ が 0 (CALC 1e-12、追跡 1e-9)、z = 0 の粒子は −V[sin(φ + ωDZ/v) − sin φ] (2e-4)。gs (TRPT、領域の外の空洞 2 台): CALC の DZ・DDP と追跡 (z = DZ と z = 0 の粒子) が全桁一致 (v2 の定義 = K64 の TRPT の動作)。Emittance[]: 軌道の位相が PHI = 0 なので、同期位相 φs の空洞とのシンクロトロンチューンの比 = 1/√cos φs (1e-3) | `test_autophi_locks_cavity_phase_to_design_orbit_arrival`、`test_emittance_autophi_uses_the_cavity_phase_at_the_orbit` | 従来どおり (NotImplementedError) |
| 2d. 領域内の RF つき MULT の Emittance[] | tmulte を bz つきで | 実装 ([多極つきの本体の半分 → 薄い空洞 → 半分] の分割)。領域内の CAVI の Emittance[] に端の RF のフリンジを足した (担当 C の実装には無く、追跡と 6e-5 チューンが違った)。bz = 0 の領域の中の CAVI・MULT は領域の外と同じ写像 | 多極なしの RF MULT = CAVI (1e-9)。bz → 0 の領域内の CAVI = 領域の外の CAVI (チューン 1e-8)。K1・K2 つきは追跡 (RFSW、放射なし) の 1 周の 6×6 行列の固有値とチューンが 3e-6 | `test_emittance_with_rf_mult_in_region` | 従来どおり |
| 3a. スピン追跡 (POL + RAD、track6 の経路): 領域内の BEND・MULT・CAVI、SOL のずれ・回転・F1 | tsol.f は領域内のすべての要素でスピンを追う (tradk の中の sprot)。pysad v1 は領域内は DRIFT・QUAD だけ | 実装。写像は高速追跡のカーネル、放射とスピンの歳差は track6 の tradk (Python)。SOL の端の回転はスピンも回す (s ← R s) | Thomas-BMT dS/dl = S × [(1+aγ)B⊥ + (1+a)B∥]/\|p\| と運動方程式の数値積分: DRIFT 5e-7、BEND 3e-6、BEND + bz 1e-5、QUAD + bz 3e-6、MULT (K0+K1) + bz 5e-5 (1 GeV、0.4 m、軸上と軸外・δ = ±5e-3、3 つの向き)。スピンと運動量の角 = aγ × 偏向角 (2e-4、偏向 0.02 と 0.3 rad)。R = 幾何の軸の回転 G1ᵀG0 (1e-12)。座標は高速追跡と同じ (1e-10) | `test_v2_C2_spin.py::test_spin_in_sol_region_matches_thomas_bmt_integration` (5)、`test_spin_precession_in_dipole_…`、`test_spin_frame_rotation_…`、`test_spin_path_follows_the_same_orbit_as_fast_tracking` (2) | 従来どおり |
| 3b. リングのスピンチューン | — | — | 平面リング (0.5 GeV、aγ = 1.135) + 長さ 1 m のソレノイド (0、1、2.5 T): 1 周のスピンの回転角 μ が cos(μ/2) = cos(π aγ) cos(φ/2)、φ = (1+a) bz L と 2e-5 | `test_ring_spin_tune_with_solenoid` (3) | — |
| 3c. UND・PHSROT のスピン追跡 | UND は磁場がなく (B37)、PHSROT は動かない (B38)。pysad は track6 の経路で NotImplementedError | UND: 刻みごとに [写像 → 放射とスピンの歳差]。縦磁場 Bs (軸外) と、周期の途中で終わる端のポテンシャルの飛び (薄い蹴り) の歳差も入る。PHSROT: 位相空間の写像 (減衰・拡散) だけで、スピンは動かさない (磁場を持たないので歳差が決まらない) | UND (1 T、1 GeV): B = rot A の中の Thomas-BMT の積分と 2e-6 (POLE = 0.25、5、5.25、中点面の外 2〜3 mm)。座標は高速追跡と丸めまで同じ | `test_spin_in_undulator_matches_thomas_bmt_integration` (4) | 従来どおり |
| N6 (新規、K64) sprot の縦磁場 bsi | `bz = bz0/pzm + bsi`: bsi = bz·al は p0 で規格化したまま z 軸の向きに足す。(1) δ≠0 で (1+a)B∥ の歳差が 1/(1+δ) 倍にならない (δ = −5e-3、0.12 rad で 6e-4 rad)、(2) 軌道が傾くと bsi·ẑ の運動量に垂直な成分が B⊥ に二重に入る、(3) 双極 + ソレノイド (tdrift) の `bsi0 − ak0x·y − ak0y·x` は T-BMT と 3.6e-3 rad 合わない | SOL 領域と UND では、運動量の向きの磁場の積分 ∫B∥dl/(Bρ) = [(bz pz + Bx px + By py)·al/pz + ψ]/(1+δ) を運動量の向きに置く | 上の T-BMT の積分 (v2 1e-7〜1e-5、K64 の式 6e-4〜4e-3) | 上と `test_k64_mode_spin_path_is_unchanged` | 再現 |
| N7 (新規、K64) tsolqur の縦磁場と放射の基準 | 領域内の四極 (RAD): bsi を最初のスライスは半分、最後のスライスは半分で上書きするので 1 スライス分足りない (1 スライスなら歳差が半分)。放射の基準の向き pxr0 が正準運動量のまま (tquad → tsolqur) で、最初のスライスに (bz/2)(y, −x) の見かけの偏向が入る | スライスごとに bz·aln/(1+δ) (と端の ±k x y/(1+δ))。放射の基準は力学的運動量 (高速追跡も) | 軸上の粒子のスピンの回転 (1+a) bz L: K64 は半分 (0.0600 / 0.1199)、v2 は 5e-15。軸外の粒子の放射損失: v2 の四極 −1.8130e-10 と MULT (K1。担当 C が放射の式と照合済み) −1.8135e-10 (3e-4)、K64 は −4.00e-10 (2.2 倍) | `…thomas_bmt_integration[quad_sol]` | 再現 |
| N8 (新規、K64) sprot の刻みの角 | 歳差の角が (1+aγ)·2 sin(θ/2) (外積 sin θ を中点の pz ≈ cos(θ/2) で割る)。刻みが粗いと θ²/24 小さい (0.26 rad の BEND を 3 刻みで 3e-4) | 厳密な θ (一様磁場の刻みで厳密)。すべての要素の追跡のスピンに効く | 平面リングの 1 周の回転角 = 2π aγ (上の 3b。K64 の式は 4e-3 rad ずれる) | `test_ring_spin_tune_with_solenoid[0.0]` | 再現 |
| 4a. RADTAPER の閉軌道の残り | tcod の収束判定は Σ(Δ/codw)² < 1e-4 (codw_x = 1e-6 m): 軌道の残りの誤差 1e-8 m。gs の既定では x = −2.19771e-7、δ = 3.5e-9 | pysad (`radtaper.solve`、判定 1e-10) の x = −2.17319317e-7 は、gs で `CODCONV = 1e-4` (判定を 1 万倍きつく) にした値 −2.17319735e-7 と 6 桁、U0 10 桁、DDP 9 桁、チューン全桁で一致: 収束した解。担当 C の「1〜2% の差」は gs の既定の判定の残り | 判定を 1e-13 にしても δ の分布が 1e-11 以内、不動点の残り 1e-12 | `test_v2_C2_optics.py::test_radtaper_closed_orbit_is_the_converged_gs_solution`、`test_radtaper_solution_is_a_fixed_point_and_tolerance_independent` | 同じ |
| 4b. RADTAPER の後の追跡 (新規、pysad) | 収束した解は閉軌道の z = 0 (trf0 = 0.038358122) | 閉軌道の z を RF の位相の基準 trf0 に移す (z → 0、trf0 → trf0 + z/v0)。担当 C の実装は同じ物理の解を z = 2.15e-5、trf0 = 0.038336637 と分けていて、z = 0 から始める追跡の位相が 2e-5 m ずれた (追跡の 1% の差の原因) | gs (CODCONV = 1e-4) の追跡 2 周: FLUC x = 4.808978e-5 (pysad 4.808984e-5)、NOFLUC・その後の NORADTAPER も 1e-5〜1e-7 | `test_radtaper_closed_orbit_is_the_converged_gs_solution` | 同じ |
| 4c. CALC の RADTAPER | 4 次元の CALC は放射を入れない (テーパーも無い: qtwiss に rtaper の式は無い)。CALC6D (qcell6d → temit) は rt = radtaper を渡すが、gs は `qcell6d---> Unstable optics` を出し、入口の軌道がテーパーなしの値 (−7.887e-5) のまま (動作しない) | CELL + CALC6D の CALC に RADTAPER を入れる (Emittance[] と同じ閉軌道・同じ係数)。4 次元の CALC は K64 と同じくテーパーなし (放射のない光学なので運動量が下がらない) | CALC6D の DX(1) = Emittance[] の閉軌道 (1e-6)、DDP = Emittance[] の DDP、「Unstable」を出さない | `test_radtaper_in_calc6d` | 従来どおり (CALC6D にテーパーなし) |
| 8. N3: CalculateOptics の TRPT の加速 | CALC と同じ結果 | CALC と同じ (TRPT かつ RFSW の加速を区間に入れる。区間の途中からでも) | gs: 終端の BX AX BY AY EX EPX NX NY が全桁、CalculateOptics − CALC = 0 | `test_calculate_optics_includes_trpt_acceleration` (2) | 同じ (pysad の不具合の修正) |
| N9 (新規、pysad) 加速区間の分散 | qmat2twiss: gr = gammab(出口)/gammab(入口) を入口の分散に掛ける | そのとおり (pysad は逆数 gammab(入口)/gammab(出口) を掛けていて、入口に分散のある加速区間の EX EPX が gs と合わなかった: 0.1109 / 0.0734) | gs と全桁 (薄い空洞・厚い空洞)。薄い空洞では x は変わらず δ = Δp/p が p0/p1 倍 → 分散は (p1/p0)^1.5 倍 (座標の規格化込み) | `test_dispersion_through_accelerating_cavity_matches_gs` (2) | 同じ (pysad の不具合の修正) |
| 10. 描画の試験 | — | `test_plot.py`: Tk がマウスの下の項目に付ける `current` タグを読み取りで外す (判定は同じ) | — | `test_plot.py` | 同じ |

担当 C の「要判断」への決定 (親の判断) の反映:

| 要判断 | 決定 | 反映 |
|---|---|---|
| B17 (NaN の返し方) | 現状のまま | 変更なし |
| DELCOR の corlength≠0 の規格化 | rms = delx が自然なら v2 で直す | 直した ([v2_changes_C2_misc.md](v2_changes_C2_misc.md))。`test_v2_C_delcor.py` の ⟨\|c\|²⟩ の期待値を 2 に更新 |
| DELCOR の CAVI の符号 | v2 の +delx | 変更なし (決定済みと記録) |
| UND の光学 | 現状 | 変更なし |
| 空間電荷の k64 モード | K64 の動作のまま | 変更なし。例 `flag_SPAC` を作り直した ([v2_changes_C2_spch.md](v2_changes_C2_spch.md)) |
| MEASURE の TRPT、`TURNS=-1` | v2 の実装でよい。v1 との違いを明記 | [v2_changes_C2_measure.md](v2_changes_C2_measure.md) §4 の表 |

## 2. 使い方・使える範囲と条件・使えない範囲

### 2.1 ずれて傾いたソレノイド (項目 1)
```
SOL SA = (L=0 BZ=2 BOUND=1 GEO=1 DX=0.002 DY=-0.001 DPX=0.003 DPY=0.002 F1=0.1) SB = (L=0 BZ=0 BOUND=1 F1=0.1);
LINE … (… SA DS QS1 DV QS2 DS SB …);
```
- **GEO=1 の SOL** (領域の始端か終端のどちらか一方) に、設計軌道の位置 `DX DY` [m] と向き `DPX DPY` [rad] をソレノイドの系で書く。
  もう一方の SOL の DX DY DZ と両端の CHI1..3 は幾何 (tsgeo) が決める (利用者が書いた値は使わない。K64 も上書きする)。
  GEO が無い (か両端にある) 領域は、SOL に書いた `DX DY DZ CHI1 CHI2 CHI3` をそのまま座標変換に使う (K64 と同じ。設計軌道が閉じるかは利用者の責任)。
- 領域の中の要素 (DRIFT・QUAD・BEND・MULT・CAVI・SOL) の軸は **ソレノイドの軸**。領域の中の Twiss・軌道 (`Twiss["DX", …]`) はソレノイドの系の値。
  領域の外では設計軌道 = 0 のまま。領域の長さ (L の和) は軸に沿った長さで、設計軌道の行路長との差は出口の SOL の DZ に入る。
- `F1` (≠0) は端の放射に使う長さ。RAD のとき、端の蹴りを長さ F1 の一様な偏向として放射させる (`RAD=1` を SOL に書くと放射しない)。
- 使える: `TrackParticles` (高速追跡、POL + RAD のスピン追跡)、CALC (CELL・INS) の Twiss と軌道、Emittance[]・EMIT、幾何 (`LINE["GX"…]`)、幾何の fit、
  MEASURE (追跡と同じ経路)。
- 使えない: k64 モード。CALC6D・RADCOD つきの CALC (SOL があると従来どおり SADError)。追跡を領域の途中から始める・途中で終える
  (`TrackParticles` の区間の端が領域の中: `ソレノイドの終端がない` のエラー。従来どおり)。
  領域より上流の幾何に軸まわりの回転 (GCHI3≠0) があるラティスで、追跡を先頭以外から始めた場合 (区間ごとに幾何を作るので、CHI3 が入らない)。
  領域内の APERT のスピン追跡の経路 (NotImplementedError)。
- 未確認: ソレノイドの軸が設計軌道と大きく (10 mrad 以上) ずれる場合、領域内の SOL (BOUND=0) にずれ・回転を書いた場合 (K64 は使わない。pysad も無視)、
  GEO のない領域に手で書いた CHI の符号の約束 (ソースの式どおりだが、gs と比べていない)。gs との数値比較は行っていない (gs は SOL + 四極 + RAD で保存値と食い違う: docs/44 §7)。

### 2.2 領域内の空洞 (項目 2)
- TRPT の加速: `TRPT; RFSW;` で、領域内の CAVI・VOLT≠0 の MULT が加速する。CALC (INS) の Twiss・軌道と `TrackParticles` に入る。
  `LINE["GAMMABETA", k]` が各位置の基準運動量。加速の後の領域の要素は、bz をその場所の運動量で規格化する。
- wake: `WakeFunction[Longitudinal|Transverse, "名前"] = {{z, w}, …}` (か旧式の `WAKE` コマンド) を領域内の CAVI・RF つき MULT に付けられる。
  wake が働くのは TRPT のとき (従来どおり)。
- AUTOPHI: `CAVI … AUTOPHI=1` (MULT も)。空洞の位相 PHI + DPHI を、**設計軌道がその空洞に着く時刻** (Twiss の DZ) を基準に測る。
  追跡は直前の CALC の DZ を使う (CALC をしていなければ CALC してから)。TRPT の CALC では設計軌道がいつも公称の加速 (VOLT·sin φ) を受ける。
  リング (Emittance[]) では軌道が位相 PHI で通る (放射損失を補うには PHI を同期位相にする)。
- RF つき MULT の Emittance[]: bz≠0 の領域の中だけ対応 (bz = 0 の領域と領域の外の RF つき MULT の Emittance[] は従来どおり未対応)。
- 使えない: k64 モード。スピン追跡の経路 (track6) での TRPT の加速 (RF は入るが基準運動量は変えない)。`V1 V20 V11 V02` (領域内の CAVI)、ANGLE≠0 の MULT (領域内)。
- AUTOPHI の意味は K64 のソースの読み (SADHelp に説明がない)。TRPT では gs と全桁一致。リング (TRPT でない) の追跡の AUTOPHI は K64 と違う (K64 は位相の差が消える: N11)。
- 未確認: 領域内の空洞の wake の横の蹴りを gs と比べていない。加速 + RAD の組み合わせ。

### 2.3 スピン追跡 (項目 3)
- `POL; RAD;` の `TrackParticles` (ビームは 9 列: 7 列目 = atan(sz/sx)、8 列目 = sy)。スピンは RAD が ON のときだけ動く (K64 と同じ)。
- v2 で通る要素: SOL 領域 (ずれ・回転・F1、領域内の DRIFT・QUAD・BEND・MULT・CAVI・SOL)、UND (磁場つきも)、PHSROT。
- 精度: 歳差は放射の刻みごとに 1 回の回転にまとめる。一様磁場の刻みは厳密、四極・多極の中は刻みの粗さ (上の表の 1e-6〜5e-5)。
  刻みを細かくするには要素の `EPS` を小さくする。
- PHSROT はスピンを動かさない。リング全体を PHSROT 1 個で表すとスピンチューンは 0 になるので、偏極の計算には使えない。
- 使えない: k64 モード (領域内は DRIFT・QUAD だけ、UND・PHSROT は NotImplementedError)。領域内の APERT。
- 限界: **Emittance[]・EMIT の POL の量 (SpinTune・Polarization など) は、SOL 領域の中は K64 の式 (N6 の bsi) のまま**。
  軸上のソレノイドでは追跡と 2e-5 で合うが、ずれて傾いたソレノイド (軌道が軸を外れる) ではスピンチューンで 2.3e-4、n0 の向きで 4e-3 rad 違う
  (`spin.tradke_spin` の `xpzb = xpz + (bsir0 + bzh·2·al)·pr²` を直すには、6 次元の微分 dxpzb も一緒に直す必要がある。未着手)。
  領域の外の Emittance[] のスピンチューンは N8 (刻みの角) の分だけ追跡 (v2) と違う (平面リングで 2.8e-5)。

### 2.4 RADTAPER (項目 4)
- 使い方は担当 C の記録 §2.5 のまま。v2 では `CELL; CALC6D; RFSW; RAD; RADCOD; RADTAPER; CALC;` でも効く。
- 閉軌道の z は 0 になり、位相の基準は `EquilibriumPosition` (trf0) に入る。
- gs と比べるときは gs 側で `CODCONV = 1e-4;` を入れる (既定の判定では軌道の残りが 1% 違う)。
- 追跡の閉軌道は Emittance[] の閉軌道と 1e-6 m 程度違う (テーパーなしでも 4e-7 m。追跡と tturne の放射の刻みの違い。K64 も同じ)。

## 3. 新しく見つけたもの (docs/21 への採番は未)

- **N6〜N8 (K64、スピン追跡)**: 上の表。SOL 領域・UND の中は v2 で直した。N8 (刻みの角) はすべての要素に効く。
- **N9 (pysad)**: 加速区間の分散の gr の向き。モードによらず直した。
- **N10 (K64)**: 空洞の中の wake の重み。v2 で直した (k64 は再現)。
- **N11 (K64、AUTOPHI)**: tcav.f・tmulti.f の `phis=merge(0.d0,merge(phic,w*trf0,autophi),trpt)` は、TRPT でないとき phis = phic にするので dphis = 0 になり、
  足したばかりの DZ·w も PHI も効かない。光学 (tcave.f: phii = phic) とも合わない。v2 は TRPT のときの形 (phis = 0) にそろえた。
  tcave.f は autophi のとき dvcacc・ddvcacc を足さないので、空洞がすべて AUTOPHI だと実効の RF の傾きが 0 になる (v2 は軌道の位相での傾きを足す)。
- **N12 (pysad、未修正)**: `Twiss["DETR"]` が CALC の後で常に 0 (`ffscmd._rows_from_tws` が detr を渡していない)。MEASURE 担当の指摘。
  結合の強いリングの MEASURE (モードが入れ替わる形) に効く。
- **N13 (K64、未修正、参考)**: tsol.f の trots の入口は z に DZ だけを足し、出口は DZ·(1 − dv) を足す。v2 は両方 (1 − dv)。
  入口の DZ≠0 は GEO が出口側の SOL にあるときだけで、差は dv·DZ (3 GeV では 1e-13 m 以下)。
- 追跡 (RAD、NOFLUC) を Emittance[] の閉軌道から始めても、1 周で x が 4e-7 m、z が 2e-6 m 動く (テーパーなし)。K64 の追跡と tturne の放射の入れ方の違いで、
  pysad の k64 の追跡は gs と 1e-12 で合っている。直していない。

## 4. 変更したファイル (C2 本体。分担の 3 件は各記録)

- `pysad/track6n_ext.py`: `trots1`・`k_sol_edge` (K_SOL_EDGE = 17)、`compile_sol_region` (v2 の端・K1=0 の四極・p0_ref)、`sol_mult_full` (加速・wake・AUTOPHI・放射のフック)、
  `step_kind` (bz の規格化・戻り値・wake)、`_compile_sol_v2` (AUTOPHI)、`compile_mult` (AUTOPHI)、`mult_full` (AUTOPHI・wake の重み)、`k_sol_quad` (放射の基準)。
- `pysad/track6n.py`: `AUTOPHI_DZ`・`autophi_dz`、`Prog(l0=)`・`autophi`、`turn` (AUTOPHI、領域内の空洞の wake と加速)、`cav_full` (wake の重み)、`TrackRun.segment`。
- `pysad/track6.py`: `_sprot` (bpar・刻みの角)、`_tsolqu_bz_ens`・`_tquad_sol` (v2)、`_tsol_region_v2`・`_trots`・`_sol_edge`・`_tradk_sol`・`_sol_drift_v2`・`_sol_bend_v2`・
  `_RadHook`・`_spin_only`・`_und_b`・`_und_v2`・`_phsrot_v2`、`_step_el` (UND・PHSROT)。
- `pysad/sgeo.py`: `region_frames`・`rot33`・`_kernel_step`、`_tsgeo1` (BEND・MULT・CAVI・K1=0 の四極)。`pysad/geo.py`: `GeoResult.regions`。
- `pysad/solmap.py` (p0・trpt・AUTOPHI)、`pysad/solemit.py` (K_SOL_EDGE、RF つき MULT、端の RF のフリンジ、bz = 0 の領域の CAVI・MULT)、
  `pysad/trptacc.py` (領域内の加速、AUTOPHI)、`pysad/lattice.py` (`_ins_twiss_acc`: gr の向き、領域の bz)、`pysad/emit6d.py` (`_z_into_trf0`、`_tcave` の AUTOPHI)、
  `pysad/emit_extra.py`・`pysad/calc6d.py` (z を trf0 に)、`pysad/match6d.py` (領域内の変数の幾何の微分)、
  `pysad/ffscmd.py` (CalculateOptics の加速、CALC6D の RADTAPER、`_z_into_trf0` の印、AUTOPHI の DZ)。
- 試験: `test_v2_C2_optics.py`・`test_v2_C2_sol.py`・`test_v2_C2_spin.py`・`test_v2_C2_solrf.py` (新規)、`test_plot.py` (current タグ)、
  `test_v2_C_measure.py`・`test_v2_C_delcor.py`・`test_v2_C_track.py` (C2 の変更に合わせて期待値 3 か所)。
- 例: `example/v2_C/C2a_sol_tilt`・`C2b_sol_rf`・`C2c_spin`・`C2d_taper_optics` (`.sad` と `.pysad.txt`)、`README.md` に追記。

## 5. 試験

| ファイル | 件数 | 内容 |
|---|---|---|
| `test_v2_C2_optics.py` | 7 | N3 (2)、加速区間の分散 (2)、RADTAPER (3) |
| `test_v2_C2_sol.py` | 11 | ずれて傾いたソレノイド (ローレンツ力 5、直進の極限 1、BEND + MULT 1)、F1 (1)、リング (2)、幾何の fit (1) |
| `test_v2_C2_spin.py` | 17 | T-BMT (領域 5、UND 4)、aγ (1)、リングのスピンチューン (3)、軌道の一致 (2)、端の回転 (1)、k64 (1) |
| `test_v2_C2_solrf.py` | 9 | 加速 (2)、CALC = 追跡 (1)、wake (2 + 1)、AUTOPHI (2)、RF MULT の Emittance (1) |
| `test_v2_C2_spch.py` | 35 | 空間電荷 (別記) |
| `test_v2_C2_measure.py` | 20 | MEASURE 系 (別記) |
| `test_v2_C2_wakecmd.py`・`test_v2_C2_delcor.py`・`test_v2_C2_script.py` | 6・7・18 | WAKE の表・DELCOR・台本 (別記) |
