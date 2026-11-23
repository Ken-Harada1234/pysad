# example/16_spin_polarization: 移動時に台本の冒頭コメントから取り除いた履歴の記述

v2_next・v2_C・v2_B の作業フォルダから機能別のフォルダへ移したとき、台本 (`.sad`・`.py`) の冒頭コメントにあった作業符号 (「v2next A1 項目 4」など) と、旧版の動作を語る文 (「旧 pysad は…」など) を取り除いた。取り除いた行を、元の行と直したあとの行の対で残す。

- `pol_fluc_tracking.py` 2 行目:
  - 元: `"""v2next A2 項目 2(b): v2 の Polarization の脱偏極率を、FLUC つきのスピン追跡 (track6、粒子ごとの Python 追跡) で測って比べる。`
  - 後: `"""v2 の Polarization の脱偏極率を、FLUC つきのスピン追跡 (track6、粒子ごとの Python 追跡) で測って比べる。`
- `pol_mpw.sad` 109 行目:
  - 元: `! ---- v2next A2 項目 2・3: MPW (逆向きの極) 入り ASTRID2 (0.58 GeV、周長 45.7 m) の RADINT と Polarization ----`
  - 後: `! ---- MPW (逆向きの極) 入り ASTRID2 (0.58 GeV、周長 45.7 m) の RADINT と Polarization ----`
- `spin_turn.sad` 1 行目:
  - 元: `! v2next A2 項目 1: example/07_spin_solenoid_special の flag_POL_RADPOL で、500 周後のスピン角が gs・k64 (0.8794 rad) と v2 (1.1516 rad) で違う原因。`
  - 後: `! example/07_spin_solenoid_special の flag_POL_RADPOL で、500 周後のスピン角が gs・k64 (0.8794 rad) と v2 (1.1516 rad) で違う原因。`
- `spin_sol_bigbz.py` 2 行目:
  - 元: `"""v2next B1 項目 6: SOL 領域内の POL (スピン)・加速 + RAD・bz ≥ 3 1/m。`
  - 後: `"""SOL 領域内の POL (スピン)・加速 + RAD・bz ≥ 3 1/m。`
- `pol_sideband.sad` 1 行目:
  - 元: `! v2next B3: シンクロトロン側帯 (スピンの位相変調) を入れた Polarization。`
  - 後: `! シンクロトロン側帯 (スピンの位相変調) を入れた Polarization。`
- `radpol_roundoff.py` 2 行目:
  - 元: `"""v2next B3 項目 2: ASTRID2 (0.58 GeV) の RADPOL 追跡で、P=1 からの緩和が ST の値の 1.58 倍になる原因 = 倍精度の丸め。`
  - 後: `"""ASTRID2 (0.58 GeV) の RADPOL 追跡で、P=1 からの緩和が ST の値の 1.58 倍になる原因 = 倍精度の丸め。`
- `pol_sideband_spin_mc.py` 2 行目:
  - 元: `"""v2next B3 項目 1: 側帯つきの脱偏極率 (spin.sideband_depolarization) を、線形の軌道 + 厳密なスピン回転の Monte Carlo で検証する。`
  - 後: `"""側帯つきの脱偏極率 (spin.sideband_depolarization) を、線形の軌道 + 厳密なスピン回転の Monte Carlo で検証する。`
- `spin_accel.py` 1 行目:
  - 元: `# v2next D1 項目 5: POL + RAD + TRPT + 加速 (track6 の経路) のスピンの歳差が、断熱加速の解析式 Φ = aγ1θ1 + aγ2θ2 になる。`
  - 後: `# POL + RAD + TRPT + 加速 (track6 の経路) のスピンの歳差が、断熱加速の解析式 Φ = aγ1θ1 + aγ2θ2 になる。`
- `pol_kick_response.py` 2 行目:
  - 元: `"""v2next E 項目 3: Polarization の脱偏極率と FLUC 追跡の差 (B3: 4.2 GeV・12 MV で追跡が式の 2.3 倍) の切り分け。`
  - 後: `"""Polarization の脱偏極率と FLUC 追跡の差 (B3: 4.2 GeV・12 MV で追跡が式の 2.3 倍) の切り分け。`
- `spin_sol_tilt.py` 2 行目:
  - 元: `"""v2next E 項目 1: ずれて傾いたソレノイド (BZ = 2 T、長さ 0.8 m の DRIFT 2 本) のあるリングの Emittance[] の POL の量 (SpinTune・n0・srot) を、追跡の 1 周のスピン行列と比べる。`
  - 後: `"""ずれて傾いたソレノイド (BZ = 2 T、長さ 0.8 m の DRIFT 2 本) のあるリングの Emittance[] の POL の量 (SpinTune・n0・srot) を、追跡の 1 周のスピン行列と比べる。`
- `spin_bend_vertical_angle.py` 2 行目:
  - 元: `"""v2next P1: 偏向磁石の中で鉛直方向の傾き py ≠ 0 の粒子のスピン回転を、厳密な Thomas-BMT の積分 (実験室系の一様磁場、出口は基準軌道に垂直な面) と、`
  - 後: `"""偏向磁石の中で鉛直方向の傾き py ≠ 0 の粒子のスピン回転を、厳密な Thomas-BMT の積分 (実験室系の一様磁場、出口は基準軌道に垂直な面) と、`
- `pol_chromatic.py` 2 行目:
  - 元: `"""v2next P1 項目 1: 4.2 GeV・12 MV の 12 セル FODO (pol_sideband.sad と同じ) の Polarization の脱偏極率に、ベータトロンのクロマティックな変調を入れる前後。`
  - 後: `"""4.2 GeV・12 MV の 12 セル FODO (pol_sideband.sad と同じ) の Polarization の脱偏極率に、ベータトロンのクロマティックな変調を入れる前後。`
- `pol_mpw_rotated.sad` 109 行目:
  - 元: `! ---- v2next P1 項目 2: 水平磁場の挿入光源 (MPW を ROTATE=90° にした ASTRID2 + MPW、0.58 GeV、45.7 m) の Polarization ----`
  - 後: `! ---- 水平磁場の挿入光源 (MPW を ROTATE=90° にした ASTRID2 + MPW、0.58 GeV、45.7 m) の Polarization ----`
- `pol_reversed_bend.sad` 1 行目:
  - 元: `! ---- v2next P1 項目 3: 反転した BEND (曲率の符号が違うシケイン) のあるリングの Polarization と PolarizationTime ----`
  - 後: `! ---- 反転した BEND (曲率の符号が違うシケイン) のあるリングの Polarization と PolarizationTime ----`
- `pol_vertical_bend.sad` 1 行目:
  - 元: `! ---- v2next P1 項目 2: 垂直の偏向 (ROTATE=90° の BEND の組) のあるリングの Polarization と PolarizationTime ----`
  - 後: `! ---- 垂直の偏向 (ROTATE=90° の BEND の組) のあるリングの Polarization と PolarizationTime ----`
- `pol_st_dk_check.py` 2 行目:
  - 元: `"""v2next P1 項目 2: v2 の Polarization の緩和率 (Sokolov-Ternov の厳密式) を、独立な DK の数値積分 (pysad/tests/_p1_dk.py) と比べる。`
  - 後: `"""v2 の Polarization の緩和率 (Sokolov-Ternov の厳密式) を、独立な DK の数値積分 (pysad/tests/_p1_dk.py) と比べる。`
- `spin_sol_quad_mult.py` 1 行目:
  - 元: `# v2next P4 項目 1: SOL 領域内の QUAD・MULT (K0 + K1) のスピン (POL + RAD) の Thomas-BMT の積分との差 (最大の成分)。`
  - 後: `# SOL 領域内の QUAD・MULT (K0 + K1) のスピン (POL + RAD) の Thomas-BMT の積分との差 (最大の成分)。`
- `spin_tracking_sol_und.sad` 先頭付近:
  - 元: `! C2c: スピン追跡`
  - 後: `! スピン追跡`
- `pol_chromatic.py` 先頭付近:
  - 元: `追跡 (FLUC つきスピン追跡、A2・B3、480 粒子)`
  - 後: `追跡 (FLUC つきスピン追跡 pol_fluc_tracking.py、480 粒子)`

## 手で書き直した記述
- 「以前の記録の…」「以前の pysad は…」「記録は docs/history/v2next_*.md §…」「(.pysad_k64.txt)」の言及など、複数行にまたがる記述は個別に書き直した (`.sad` 側に残るのは現在の動作の説明だけ)。

## 削除した k64 モードの出力
- 旧 `v2_next` の `pol_sideband`・`pol_mpw_rotated`・`pol_reversed_bend`・`pol_vertical_bend` の `*.pysad_k64.txt` と `spin_sol_quad_mult.out_k64.txt` を削除した (k64 モードは `PYSAD_COMPAT=k64` で再現できる)。

## 仕上げで直した台本の記述 (旧パス・作業符号の点検)
- `spin_turn.sad` 2 行目: 元 `原因: K64 の sprot の刻みの角 (1+aγ)·2 sin(θ/2) (C2 の N8)。` → 後 `原因: K64 の sprot の刻みの角 (1+aγ)·2 sin(θ/2)。`
- `pol_kick_response.py` 2〜3 行目: 元 `差 (B3: 4.2 GeV・12 MV で追跡が式の 2.3 倍)`・`線形の軌道モデル (M1: B3 の MC と同じ)` → 後 `差 (4.2 GeV・12 MV で追跡が式の 2.3 倍)`・`(M1: pol_sideband_spin_mc.py と同じ)`
- `spin_sol_quad_mult.py` 3 行目: 元 `v2 はどの bz でも 1e-6 以下 (P4 の前は QUAD 9e-5、MULT 1.6e-4)。` → 後 `v2 はどの bz でも 1e-6 以下。` (取り除いた「P4 の前」の値: QUAD 9e-5、MULT 1.6e-4)
