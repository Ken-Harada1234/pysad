# example/16_spin_polarization の gs との比較 (README の表から移した記述)

`example/16_spin_polarization/README.md` の表の「gs と pysad」の列にあった、gs (K64) の出力との比較の記述を、元の表のまま移した。README の本文は、結果の数字の説明を主役にし、gs との比較は「検証」の 1〜3 行にした。

| 名前 | 内容 | gs と pysad |
|---|---|---|
| [pol_sideband](pol_sideband.sad) | 4.2 GeV・RF 12 MV の 12 セル FODO (偏向磁石 15 度 × 12、四極の縦ずれ ±3 mm、RF 調和数 64) の Polarization。シンクロトロン側帯 (スピンの位相変調) の項 | △ (gs・k64 は P = 1.27e-4 (spdepol が脱偏極率を 1000 倍大きく出す)、v2 は P = 0.119 (側帯の項込み。線形の式だけなら 0.205)。偏極時間は gs 23.06 s、v2 23.08 s (緩和率の定義の違い)) |
| [pol_sideband_spin_mc](pol_sideband_spin_mc.py) | 側帯の式の検証。線形の軌道 + 厳密なスピン回転の Monte Carlo。`py -3.12 example/16_spin_polarization/pol_sideband_spin_mc.py 4.2 12 3000 30000 11 out.json` (3000 粒子 × 30000 周で約 8 分) | gs は使わない |
| [pol_fluc_tracking](pol_fluc_tracking.py) | v2 の Polarization の脱偏極率を、FLUC つきのスピン追跡 (track6、粒子ごとの Python 追跡) で測って比べる。`E[GeV] V[MV] 粒子数 乱数の種 n1 ndamp` (例: `4.2 12 120 11 1500 4500`。1 粒子 1 周約 4.5 ms)。結果は `dep_<E>_<種>.json` | gs は使わない |
| [pol_chromatic](pol_chromatic.py) | ベータトロンのクロマティックな変調を入れる前後の脱偏極率 (`pol_sideband.sad` と同じ 4.2 GeV・12 MV の FODO、約 1 分)。測定の保存データは [`pol_chromatic_response_samples_f1.json`](pol_chromatic_response_samples_f1.json)・[`f2_models`](pol_chromatic_response_samples_f2_models.json)。結果は `pol_chromatic.out.txt` | FLUC つき追跡の (10.4 ± 1.6)e-8 /周 (480 粒子) と比べる。gs は使わない |
| [pol_kick_response](pol_kick_response.py) | 脱偏極率と FLUC 追跡の差 (4.2 GeV・12 MV で追跡が式の 2.3 倍) の切り分け。δ の蹴りへの応答 (1 試料 1〜2 分)。結果は `pol_kick_response.out.txt` | gs は使わない |
| [pol_mpw](pol_mpw.sad) | MPW (逆向きの極) 入り ASTRID2 (0.58 GeV、周長 45.7 m) の RADINT と Polarization | △ (K64 は RADINT が K0 の極を無視 (U0 6.17 keV)、Polarization は NaN、PolarizationTime は符号つきの率 (28842 s)。v2 は U0 7.81 keV、P = 0.7019、τ = 21915 s で独立な式 (P = 0.7016、τ = 21886 s) と合う) |
| [pol_mpw_rotated](pol_mpw_rotated.sad) | MPW を `ROTATE = 90°` にした水平磁場の挿入光源 (ASTRID2 + MPW) の Polarization | △ (K64 は緩和率の Σ r (b·n) が MPW の放射を数えず τ が実際より短い。v2 は独立な DK の数値積分と τ が 2e-3、P が 4e-3 で合う (0.7007 対 0.7018、τ 21914 s)) |
| [pol_vertical_bend](pol_vertical_bend.sad) | 垂直の偏向 (`ROTATE=90°` の BEND の組、軌道を鉛直に ±0.1 rad) のある 0.5 GeV の 12 セル FODO の Polarization と PolarizationTime | △ (gs・k64 は垂直の偏向の放射を数えず Polarization が NaN。v2 は Sokolov-Ternov の厳密式で、独立な DK の数値積分と τ が 1e-5、P が 1e-3) |
| [pol_reversed_bend](pol_reversed_bend.sad) | 反転した BEND (曲げ角 +0.05・−0.10・+0.05 rad のシケイン) のある 3 GeV の 12 セル FODO の Polarization と PolarizationTime | △ (K64 は Polarization が NaN、PolarizationTime は符号つきの率。v2 は P・τ が独立な式 (P = 0.66825) と 1e-3) |
| [pol_st_dk_check](pol_st_dk_check.py) | v2 の緩和率 (Sokolov-Ternov の厳密式) を独立な DK の数値積分 (`pysad/tests/_p1_dk.py`) と比べる。0.58 GeV の ASTRID2 + MPW の MPW を `ROTATE` = 0.3 rad / 90° にしたもの (約 6 分)。結果は `pol_st_dk_check.out.txt` | gs は使わない |
| [spin_turn](spin_turn.sad) | `07_spin_solenoid_special` の `flag_POL_RADPOL` で、500 周後のスピン角が gs・k64 (0.8794 rad) と v2 (1.1516 rad) で違う原因。K64 の sprot の刻みの角 (1+aγ)·2 sin(θ/2) | △ (放射の分割 EPS を細かくすると k64 と v2 が一致する: 1 周の差 5.6e-4 rad → 8e-8 rad) |
| [radpol_roundoff](radpol_roundoff.py) | ASTRID2 (0.58 GeV) の RADPOL 追跡で、P = 1 からの緩和が Sokolov-Ternov の値の 1.58 倍速い原因 (倍精度の丸め)。出力 (`radpol_roundoff.pysad.txt`): P0 = 1、0.5、0.25 から始めた 1 周あたりの dP と式の比 (1.578、0.934、0.979) | P0 = 1 だけ増分が倍精度の ulp 未満。gs は使わない |
| [spin_tracking_sol_und](spin_tracking_sol_und.sad) | スピン追跡 (POL + RAD) の SOL 領域・UND。平面リング + ソレノイドの 1 周のスピン回転角 (cos(μ/2) = cos(π aγ) cos(φ/2))、領域内の BEND のスピンと運動量の角 (aγ × 偏向角)、UND (5 周期、1 T) でスピンが元に戻る | v2 のみ (K64 に無い機能) |
| [spin_sol_bigbz](spin_sol_bigbz.py) | SOL 領域内の DRIFT のスピン (bz = 4.5・9 1/m)・POL の SpinTune・加速 + RAD。Thomas-BMT・解析式との差の表。結果は `spin_sol_bigbz.out.txt` | gs は使わない |
| [spin_sol_tilt](spin_sol_tilt.py) | ずれて傾いたソレノイド (BZ = 2 T、長さ 0.8 m の DRIFT 2 本) のあるリングの Emittance[] の POL の量 (SpinTune・n0・srot) を、追跡の 1 周のスピン行列と比べる (約 1 分)。結果は `spin_sol_tilt.out.txt` | gs は使わない |
| [spin_sol_quad_mult](spin_sol_quad_mult.py) | SOL 領域内の QUAD・MULT (K0 + K1) のスピン (POL + RAD) の Thomas-BMT の積分との差 (最大の成分)。0.4 m の要素、軸上と軸を外れた粒子、3 つのスピンの向き。`py -3.12 example/16_spin_polarization/spin_sol_quad_mult.py [v2|k64]` | v2 は 1e-6 以下 (`spin_sol_quad_mult.out.txt`)。gs は使わない |
| [spin_bend_vertical_angle](spin_bend_vertical_angle.py) | 偏向磁石の中で鉛直方向の傾き py ≠ 0 の粒子のスピン回転を、厳密な Thomas-BMT の積分と、Emittance[] (tturne) の srot・追跡 (track6) で比べる。扇形・K1 = 0・L = 1.2 m・15°・0.5 GeV (数十秒)。結果は `spin_bend_vertical_angle.out.txt` | gs は使わない |
| [spin_accel](spin_accel.py) | POL + RAD + TRPT + 加速 (track6 の経路) のスピンの歳差が、断熱加速の解析式 Φ = aγ1θ1 + aγ2θ2 になる。1 GeV、扇形の偏向磁石の間に加速する空洞 (VOLT = 0・0.2・0.5 GV)。結果は `spin_accel.out.txt` | v2 のみ (k64 は SADError) |
