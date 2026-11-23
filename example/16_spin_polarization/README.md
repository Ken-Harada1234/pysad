# 16_spin_polarization (スピン追跡・偏極度・緩和率)

スピンの追跡 (`POL`・`RAD`)、Polarization・PolarizationTime (緩和率)、シンクロトロン側帯、クロマティックな変調、ソレノイド領域内のスピン、加速中のスピンの歳差。
ソレノイド領域そのものの光学・光子は `07_spin_solenoid_special` にある。ここは偏極とスピンの量を、独立な式・厳密な積分・追跡と比べる例を集めた。

この文書で分かること: 各例が何を計算し、結果の数字 (偏極度 P、偏極時間 τ、スピン角など) が何を意味するか。前提: 加速器物理と SAD の基本 (`POL`・`RAD`・`FLUC` は SAD のフラグ)。読む順: 下の「用語」→ 表 → 「検証」。

## 用語
- **Polarization** (P): 電子ビームが放射で自然に偏極する平衡の偏極度 (Sokolov-Ternov の値に、スピンと軌道の結合による脱偏極を入れたもの)。`Polarization` が返す。理想的な平らなリングで 0.92 (≈ 8/(5√3)) が上限の目安。
- **PolarizationTime** (τ [s]): 偏極が平衡に近づく時定数。`PolarizationTime[s]`。
- **脱偏極率**: 放射の量子ゆらぎがスピンの向きを乱す率 [1/周]。P を下げる。
- **側帯**: シンクロトロン振動 (周波数 νs) でスピン位相が変調されるために出る、スピンチューン ± νs の共鳴線。

## 走らせ方
`.sad` の台本は 3 ファイル (`NAME.sad`、`NAME.gs.txt` (gs の画面表示。ない例もある)、`NAME.pysad.txt` (pysad v2 の画面表示))。
`.py` の例は、リポジトリ直下から `py -3.12 example/16_spin_polarization/NAME.py` で走らせ、出力は `NAME.out.txt` (一部は引数つき。台本の先頭の「使い方」を参照)。
台本は `py -3.12 example/tools/run_example.py example/16_spin_polarization/NAME.sad` (gs と pysad の両方) か、`PYTHONPATH=.` で `pysad.ffscmd.run_sad_string` に台本を渡す。既定は v2。gs と同じ動作で走らせるときは環境変数 `PYSAD_COMPAT=k64` (k64 モード)。

## 例の一覧と結果の要点
表は、例ごとに「何の例か」と、結果の数字が何を意味するかを並べる。

| 名前 | 内容 | 結果の要点 |
|---|---|---|
| [pol_sideband](pol_sideband.sad) | 4.2 GeV・RF 12 MV の 12 セル FODO (偏向磁石 15 度 × 12、四極の縦ずれ ±3 mm、RF 調和数 64) の Polarization。シンクロトロン側帯 (スピンの位相変調) の項 | v2 は P = 0.119 (側帯の項込み。線形の式だけなら 0.205)、偏極時間 τ = 23.08 s。縦ずれがあるので脱偏極が効き、理想値 0.92 より大きく下がる |
| [pol_sideband_spin_mc](pol_sideband_spin_mc.py) | 側帯の式の検証。線形の軌道 + 厳密なスピン回転の Monte Carlo。`py -3.12 example/16_spin_polarization/pol_sideband_spin_mc.py 4.2 12 3000 30000 11 out.json` (3000 粒子 × 30000 周で約 8 分) | 脱偏極率は 5.0e-8 /周 (側帯の式 4.5e-8、線形の式 2.3e-8)。gs は使わない |
| [pol_fluc_tracking](pol_fluc_tracking.py) | v2 の Polarization の脱偏極率を、FLUC つきのスピン追跡 (track6、粒子ごとの Python 追跡) で測って比べる。`E[GeV] V[MV] 粒子数 乱数の種 n1 ndamp` (例: `4.2 12 120 11 1500 4500`。1 粒子 1 周約 4.5 ms)。結果は `dep_<E>_<種>.json` | 4.2 GeV・12 MV で脱偏極率 (10.4 ± 1.6)e-8 /周 (480 粒子)。gs は使わない |
| [pol_chromatic](pol_chromatic.py) | ベータトロンのクロマティックな変調を入れる前後の脱偏極率 (`pol_sideband.sad` と同じ 4.2 GeV・12 MV の FODO、約 1 分)。測定の保存データは [`pol_chromatic_response_samples_f1.json`](pol_chromatic_response_samples_f1.json)・[`f2_models`](pol_chromatic_response_samples_f2_models.json)。結果は `pol_chromatic.out.txt` | FLUC つき追跡の (10.4 ± 1.6)e-8 /周 (480 粒子) と比べる。gs は使わない |
| [pol_kick_response](pol_kick_response.py) | 脱偏極率と FLUC 追跡の差 (4.2 GeV・12 MV で追跡が式の 2.3 倍) の切り分け。δ の蹴りへの応答 (1 試料 1〜2 分)。結果は `pol_kick_response.out.txt` | 光子 1 個の δ の蹴りに対するスピンの不変軸の応答を、追跡と線形の軌道モデルで測って比を出す。差の原因の絞り込みで、完全な説明には至っていない (「未確認・限界」) |
| [pol_mpw](pol_mpw.sad) | MPW (逆向きの極) 入り ASTRID2 (0.58 GeV、周長 45.7 m) の RADINT と Polarization | v2: U0 = 7.81 keV、P = 0.7019、τ = 21915 s。独立な式の P = 0.7016、τ = 21886 s と合う。逆向きの極があるので、放射の積分に K0 の極を入れる必要がある |
| [pol_mpw_rotated](pol_mpw_rotated.sad) | MPW を `ROTATE = 90°` にした水平磁場の挿入光源 (ASTRID2 + MPW) の Polarization | 磁場が水平になると b·n = 0 で、軌道は鉛直に波打つ。v2: P = 0.7007 対 独立な DK の数値積分 0.7018 (4e-3)、τ 21914 s (2e-3) |
| [pol_vertical_bend](pol_vertical_bend.sad) | 垂直の偏向 (`ROTATE=90°` の BEND の組、軌道を鉛直に ±0.1 rad) のある 0.5 GeV の 12 セル FODO の Polarization と PolarizationTime | v2 は Sokolov-Ternov の厳密式で、独立な DK の数値積分と τ が 1e-5、P が 1e-3 |
| [pol_reversed_bend](pol_reversed_bend.sad) | 反転した BEND (曲げ角 +0.05・−0.10・+0.05 rad のシケイン) のある 3 GeV の 12 セル FODO の Polarization と PolarizationTime | v2: P = 0.66825 (独立な式と 1e-3)。反転した曲げでは偏極の向きが逆の寄与が入り、P が理想値より下がる |
| [pol_st_dk_check](pol_st_dk_check.py) | v2 の緩和率 (Sokolov-Ternov の厳密式) を独立な DK の数値積分 (`pysad/tests/script_tests/_p1_dk.py`) と比べる。0.58 GeV の ASTRID2 + MPW の MPW を `ROTATE` = 0.3 rad / 90° にしたもの (約 6 分)。結果は `pol_st_dk_check.out.txt` | 緩和率が数値積分と合う。gs は使わない |
| [spin_turn](spin_turn.sad) | `07_spin_solenoid_special` の `flag_POL_RADPOL` で、500 周後のスピン角が gs・k64 (0.8794 rad) と v2 (1.1516 rad) で違う原因。K64 の sprot の刻みの角 (1+aγ)·2 sin(θ/2) | 放射の分割 EPS を細かくすると k64 と v2 が一致する (1 周の差 5.6e-4 rad → 8e-8 rad)。差は放射の分割の粗さに由来する |
| [radpol_roundoff](radpol_roundoff.py) | ASTRID2 (0.58 GeV) の RADPOL 追跡で、P = 1 からの緩和が Sokolov-Ternov の値の 1.58 倍速い原因 (倍精度の丸め)。出力 (`radpol_roundoff.pysad.txt`): P0 = 1、0.5、0.25 から始めた 1 周あたりの dP と式の比 (1.578、0.934、0.979) | P0 = 1 だけ増分が倍精度の ulp 未満で、見かけの緩和が速くなる。gs は使わない |
| [spin_tracking_sol_und](spin_tracking_sol_und.sad) | スピン追跡 (POL + RAD) の SOL 領域・UND。平面リング + ソレノイドの 1 周のスピン回転角 (cos(μ/2) = cos(π aγ) cos(φ/2))、領域内の BEND のスピンと運動量の角 (aγ × 偏向角)、UND (5 周期、1 T) でスピンが元に戻る | v2 のみ (K64 に無い機能) |
| [spin_sol_bigbz](spin_sol_bigbz.py) | SOL 領域内の DRIFT のスピン (bz = 4.5・9 1/m)・POL の SpinTune・加速 + RAD。Thomas-BMT・解析式との差の表。結果は `spin_sol_bigbz.out.txt` | 強いソレノイド磁場でも Thomas-BMT の積分・解析式と合う。gs は使わない |
| [spin_sol_tilt](spin_sol_tilt.py) | ずれて傾いたソレノイド (BZ = 2 T、長さ 0.8 m の DRIFT 2 本) のあるリングの Emittance[] の POL の量 (SpinTune・n0・srot) を、追跡の 1 周のスピン行列と比べる (約 1 分)。結果は `spin_sol_tilt.out.txt` | Emittance[] の POL の量が追跡の 1 周のスピン行列と合う。gs は使わない |
| [spin_sol_quad_mult](spin_sol_quad_mult.py) | SOL 領域内の QUAD・MULT (K0 + K1) のスピン (POL + RAD) の Thomas-BMT の積分との差 (最大の成分)。0.4 m の要素、軸上と軸を外れた粒子、3 つのスピンの向き。`py -3.12 example/16_spin_polarization/spin_sol_quad_mult.py [v2|k64]` | v2 は差が 1e-6 以下 (`spin_sol_quad_mult.out.txt`)。gs は使わない |
| [spin_bend_vertical_angle](spin_bend_vertical_angle.py) | 偏向磁石の中で鉛直方向の傾き py ≠ 0 の粒子のスピン回転を、厳密な Thomas-BMT の積分と、Emittance[] (tturne) の srot・追跡 (track6) で比べる。扇形・K1 = 0・L = 1.2 m・15°・0.5 GeV (数十秒)。結果は `spin_bend_vertical_angle.out.txt` | 厳密な積分と srot・追跡の差の表。gs は使わない |
| [spin_accel](spin_accel.py) | POL + RAD + TRPT + 加速 (track6 の経路) のスピンの歳差が、断熱加速の解析式 Φ = aγ1θ1 + aγ2θ2 になる。1 GeV、扇形の偏向磁石の間に加速する空洞 (VOLT = 0・0.2・0.5 GV)。結果は `spin_accel.out.txt` | v2 のみ (k64 は SADError) |

## 検証
- 緩和率は Sokolov-Ternov の厳密式 Σ r (1 − 2/9 (n·v)²) で、独立な DK の数値積分 (`pol_st_dk_check`)・独立な式 (`pol_reversed_bend`・`pol_mpw`) と比べた。`PYSAD_COMPAT=k64` では gs と同じ結果が出る (gs の Polarization は、垂直の偏向・反転した BEND・逆向きの極のあるリングで NaN になるか符号つきの率しか出さないので、これらの例の基準にしない)。gs との比較の詳細は `docs/history/example_16_spin_polarization_validation.md`。
- 側帯の項 (`pol_sideband`) は、FLUC つきスピン追跡 (4.2 GeV・12 MV で (10.4 ± 1.6)e-8 /周、480 粒子) と、線形の軌道 + 厳密なスピン回転の MC (5.0e-8) を基準にした。式の値は、線形の式 2.3e-8、側帯の式 4.5e-8。
- 試験は `pysad/tests/script_tests/test_v2_B3_pol.py`・`test_v2_E_pol_kick.py`・`test_v2_P1_pol_chrom.py`・`test_v2_E_spinsol.py`・`test_v2_P4_sol_spin.py` ほか。

## 未確認・限界
- 側帯の項を入れた v2 の脱偏極率は追跡 (FLUC つき) より小さいまま残る (4.2 GeV・12 MV で追跡が式の約 2.3 倍。切り分けは `pol_chromatic`・`pol_kick_response`)。
- 例のエネルギーは 0.5〜4.2 GeV、ラティスは小さな FODO と ASTRID2 (+ MPW)。実機のラティスでは確かめていない。
- `pol_sideband_spin_mc`・`pol_fluc_tracking` は統計誤差が数 % (粒子数と周数で決まる)。
