# v2 の続き (担当 B3): 偏極のシンクロトロン側帯、RADPOL の緩和、光子 (tphrec)、幾何の fit の解析微分、NAFF、MPW の偏極

作成: 2026-11-18。対象: `pysad/spin.py` (`sideband_depolarization`・`_spin_frames`)、`track6.py` (光子の記録)、`sadrandom.py` (`tdusrnpl`)、`ffscmd.py` (PHOTONS の経路)、
`match6d.py` (`tdgeo_term` の四極の DX・DY、`sol_dcod`)、`sgeo.py`・`track6n_ext.py`・`solemit.py` (SOL 領域内の四極のずれ)、新規 `pysad/naff.py`。
試験は `pysad/tests/test_v2_B3_*.py` (pol 7、radpol 3、photons 6、naff 12、geofit 4 の計 32 件。既存の `test_v2_A_spin.py` の 3 か所の数値を側帯の分だけ更新した: §2.6)。
走る例は `example/v2_next/B3_*` (§9)。互換モード: v2 の修正は `compat.v2()` で分け、k64 は K64 (gs) のまま。numba のカーネル `k_sol_quad` を変えたのでキャッシュ
(`pysad/__pycache__/track6n_ext.*.nbi・.nbc`) を消した。

## 使った条件 (リングの呼び方)

- 小さなリング: 12 セルの平面 FODO (周長 45.6 m、偏向磁石 L=1.2 m・15 度 × 12 = 180 度、四極の縦ずれ ±3 mm、RF 調和数 64)。`test_v2_A_spin.py` の `_lat` と同じ。
  K64 のモデルは周回の曲げが 180 度なので、スピンチューンは aγ/2 (4.2 GeV で aγ = 9.53 に対し端数 −0.236)。
- 24 セルの平面 FODO (3 GeV、周長 105.6 m、RF 4 MV)。大きなリング A (`pysad/tests/kek-ls_v03_12.lat`、3 GeV、周長 603 m、3563 要素)、
  SP8 (大きなリング B、8 GeV、周長 1436 m。SPring-8 型、デザインレポートより再現。数値だけ)、PF リング (`example/pf/lattice/v9_16_asym.sad`、2.5 GeV、187 m)、
  ASTRID2 (`example/astrid/lattice/ASTRID2_v01_01.sad`、0.58 GeV、45.7 m) とその MPW 入り (`ASTRID2_v01_02_id.sad`)。

## 1. 一覧

| ID | K64 の動作 | v2 の動作 | 根拠 (式・数値) | 試験 | k64 モードの扱い |
|---|---|---|---|---|---|
| SB (項目 1) | `spdepol` の側帯の項は大きさが合わない (§2.5): 4.2 GeV・12 MV の例で P = 1.27e-4、v2 (A の線形の式) は P = 0.205。追跡の脱偏極率より 3 桁大きい | 平衡偏極度の脱偏極率に、スピンの歳差の位相の変調 (シンクロトロン振動の側帯) を厳密に取った項を入れる (`spin.sideband_depolarization`)。`Polarization`・`depol_rate`・`taudep` が変わる。`depol_rate_linear` に線形の値を残す | 独立な導出 (§2.1)。(a) 変調 OFF で線形の式に 1e-5、(b) 振幅を 1/100 にした極限で線形の式に 0.6%、(c) 線形の軌道 + 厳密なスピン回転の MC (3000 粒子 × 30000 周): 4.2 GeV・12 MV の実際の振幅で MC 5.01e-8 /周、式 4.46e-8 (−11%)、線形 2.32e-8 (−54%)。RF 6 MV で MC 3.83e-8、式 3.50e-8 (−9%)。(d) FLUC つきのスピン追跡: §2.3 | `test_v2_B3_pol.py` (7 件) | 従来どおり (K64 の `spdepol`・NaN) |
| RP (項目 2) | A2: 追跡 (RADPOL) の P=1 からの緩和が (1 − PST)/τ の 1.58 倍速い | **変更なし。K64 の誤りではない**: 追跡の式は Sokolov–Ternov と厳密に一致 (3 GeV FODO で 1e-4)。ASTRID2 の 1.58 倍は、1 刻みの緩和量 6.4e-17 が倍精度の ulp (1.1e-16) 以下で、P≈1 のとき丸めで 1.73 倍になるため (§3) | 式の導出 (§3.1)、3 GeV の追跡 (1e-4)、ASTRID2 の P0 = 1、0.5、0.25 からの緩和 (比 1.578、0.934、0.979 = 丸めの量子化の予想どおり) | `test_v2_B3_radpol.py` (3 件) | 同じ |
| PH (項目 3) | BEAMBEAM 以外の光子は pysad が記録しない (BEND があると SADError) | BEND (K1=0・K1≠0・ROTATE・DX・DY・DROTATE) と QUAD の光子を `PhotonList` に記録 (tdusrnpl・radangle1・tphrec)。**K64 の偏光 xi1 は ±1 を超える誤り (1e3〜1e6) を v2 は直す** | gs と光子 358 個が全列一致 (k64 モード。E 列だけ 3e-8)。解析式 (v2): 光子数 1e-3 (統計 4e-3)、<u>/u_c = 8/(15√3)、<u²> = 11/27、スペクトルの KS 検定、位置が円弧の上 | `test_v2_B3_photons.py` (6 件) | gs と同じ (xi1 の誤りも) |
| GF (項目 4) | tdgeo は qdtrans (eps=1e-6 の中心差分) と tftmatu。pysad は領域内の四極のずれを追跡・CALC・Emittance[]・tsgeo で扱えなかった (NotImplementedError) | 領域内の四極の DX・DY を追跡のカーネル・solemit・sgeo に入れた (gs と同じ)。幾何の fit の DX・DY の微分は、(I − J)∂t/∂DX の解析式 (`match6d.sol_dcod`) と tdgeo の式 | 領域内の DX・DY・ROTATE つきの四極の幾何が gs と全桁一致。幾何の fit (QF の DX、QD の DY) が gs と同じ 2 回の Newton で収束し、変数が 1e-12 以内。解析微分は幾何の有限差分と 2e-3 | `test_v2_B3_geofit.py` (4 件) | 領域内の四極のずれは k64 でも動く (gs と同じ。SOL 自体のずれは k64 で未対応のまま) |
| NF (項目 5) | MAIN の TRACK の Fourier (wjfit・wfres・wfbin・wiord・wjdraw、FOURIER フラグ) は FFS から使えない。gs の MAIN の TRACK は SIGSEGV (`tinip`)。pysad の maintrack は FOURIER を無視し、wtunem (FFT の放物線補間) だけ | `pysad/naff.py`: NAFF (Hann 窓、Brent 法、Gram-Schmidt)。`naff`・`naff_real`・`tune_xp`・`tunes_from_coords` | 既知の周波数の信号 (1e-11)、小さなリングの追跡のチューンが CALC と 1e-9 (A 6e-10、PF 1.4e-10、ASTRID2 5e-11)。wtunem は 1e-4 | `test_v2_B3_naff.py` (12 件) | 新しい関数 (ゲートなし) |
| MP (項目 6) | A2 の P1: 逆向きの極の P・τ | 変更なし (P1 のまま)。側帯の項を入れた後の確認: ASTRID2 + MPW で P = 0.7019、τ = 21915 s が独立な式 (P = 0.70163、τ = 21886 s) と 2e-3 | A2 の値の再確認ではなく、側帯の項を入れたコードで再計算 (平面のリングなので脱偏極率は 0 のまま) | `test_v2_B3_pol.py::test_astrid2_mpw_polarization_after_sideband_v2` | NaN (B55) |

## 2. 項目 1: シンクロトロン側帯の項 (Polarization の脱偏極率)

### 2.1 導出 (SLIM の 1 次の摂動を、縦の振動数と減衰を入れて展開)

- 要素 i の出口の軌道 u (6 次元) に対し、K64 のスピンの回転は s' = (I − [D_i u]×) R_i s (D_i = Emittance[] の srot、`spin_orbit` の元)。
  閉軌道のスピン軸 n_i に垂直な基底 (m_i, l_i) を閉軌道のスピンの歳差 R_i で運び、α = s·(m + i l) とすると、α' = e^{i p_i·u} α + g_i·u、
  p_i = −n_i·D_i (n のまわりの回転 = 位相の変調、主にエネルギー δ: スピンチューンが aγ(1+δ))、g_i = i (m_i + i l_i)·D_i (n に垂直な回転 = 結合)。
- 1 周では α_{n+1} = e^{iμ}(e^{iP·u_n} α_n + G·u_n) (P = Σ p_i Φ_i、G = Σ g_i Φ_i、Φ_i は入口から要素 i への軌道の写像、μ = 閉軌道のスピンの 1 周の歳差角)。
- 放射の跳びが無いときに保存される量 α − α*(u) の不変な場 (未来の軌道についての和。減衰があるので収束する):
  α*(u) = −Σ_{k≥0} e^{−iμk} e^{−iS_k} G·M^k u、S_k = Σ_{t=0}^{k} P·M^t u (M = 減衰つきの 1 周写像)。
  変調 P を切ると α* = −G(1 − e^{−iμ}M)^{-1}u で、線形の SLIM の N と一致する。
- 放射の跳び ξ (共分散 B_i) で偏極は ½|α*(u+ξ) − α*(u)|² だけ減るので、脱偏極率 = ½ Σ_i <G_i B_i G_i^†> (G_i = 要素 i の出口での ∂α*/∂u、u はビーム行列のガウス分布で平均)。
  勾配は ∂α*/∂u = −Σ_k e^{−iμk} e^{−iS_k} [G M^k − i (G M^k u) Q_k]、Q_k = Σ_{t≤k} P·M^t。要素ごとの勾配は鎖
  G_{i+1} = [e^{i p_i u_i} G_i + i α*_i e^{i p_i u_i} p_i + g_i] T_i^{-1}、α*_{i+1} = e^{i p_i u_i} α*_i + g_i u_i で運ぶ。
- 実装 (`spin.sideband_depolarization`): M の固有分解で S_k・M^k u を k の閉じた式にし、ガウス分布の平均は Sobol 列 (1024 点、`seed` 固定) の準モンテカルロ。
  `Emittance[]` では v2・RFSW ON のとき `slim_depolarization(..., beam_matrix=…, sideband=True)` が使う。`spin.SIDEBAND = False` で線形の式だけに戻る。
  戻り: `depol_rate` (側帯込み)、`depol_rate_linear`、`spin_spread_sb` (入口での n の広がり ½<|α*|²>。`spin_spread` は線形のまま)。
  縦の減衰が無い (RF OFF) とき・不安定なリングは側帯の項を入れずに線形の式を使い、`sideband_error` に理由を残す。
- 変調の大きさ: 12 セル FODO の 4.2 GeV で、変調指数 aγσδ/νs ≈ 0.5 (aγ は周回 180 度なので 4.76、σδ = 1.84e-3、νs = 0.030)。
  スピンチューンの端数 −0.236 が νx − 2νs = 0.0547 − 0.0599 に近く、2 次の側帯の共鳴が効く。

### 2.2 検証 1: 線形の軌道 + 厳密なスピン回転の Monte Carlo (導出だけを調べる)

`example/v2_next/B3_spin_mc.py`: 要素ごとに u ← T_i u + (跳び、共分散 B_i)、s ← exp(−[D_i u]×) R_i s (D は Emittance[] の要素ごとの srot)。1 − <S·n0> の傾きが 1 周あたりの脱偏極率
(3000 粒子、30000 周、4000 周以降の直線。傾きの標準誤差は最小二乗で 0.5%、粒子の統計とあわせて ~2%)。

| 条件 | 線形の式 | 側帯の式 (1024 点) | MC (傾き) | 式/MC |
|---|---|---|---|---|
| 4.2 GeV・12 MV (νs = 0.030、実際の振幅) | 2.318e-8 | 4.456e-8 | (5.01 ± 0.10)e-8 | 0.89 |
| 同、跳びと分布を 1/100 倍 (率/f) | 2.318e-8 | 2.331e-8 | 2.386e-8 | 0.98 |
| 4.2 GeV・6 MV (νs = 0.020) | 2.330e-8 | 3.495e-8 | (3.83 ± 0.08)e-8 | 0.91 |
| 4.2 GeV・24 MV (νs = 0.043) | 2.297e-8 | 4.933e-8 | (5.50 ± 0.11)e-8 | 0.90 |

- **線形の式は 0.01 倍の振幅で MC と 3% で合い (導出の線形の極限)、実際の振幅では 2.2 倍過小評価**。側帯の式は MC と −9〜−11% で合う (RF 6・12・24 MV の 3 点で 0.91、0.89、0.90)。
- 残る −10% の原因 (確かめていない推定): 結合項 b(u) の u の 2 次の項 (1 周の途中の位相の変調 × 結合 = 相対 ~aγδθ ≈ 0.1) を入れていない。
- ASTRID2 (0.58 GeV、四極の縦ずれ rms 1 mm、60 台): 線形 1.317e-12、側帯 1.432e-12 (+8.8%)。MC (2000 粒子 × 36000 周) の傾き (1.7 ± 0.3)e-12 は統計誤差で両方と合う (率が小さく、本数を増やさないと区別できない)。

### 2.3 検証 2: FLUC つきのスピン追跡 (非線形の軌道を含む)

方法は A2 の §3.2 と同じ (`example/v2_next/A2_pol_fluc_tracking.py`: 高速追跡で 8000 周かけて平衡にし、n(u) = n0 + N0 u でスピンを置き、FLUC つきスピン追跡 1500 周と 3000 周の
1 − <S·n0> の差から率。FLUC を切って 4500 周で軌道を減衰させる。1 粒子 1 周 4.5 ms)。4.2 GeV・12 MV:

| 実行 (各 120 粒子、乱数の種) | D(1500) | D(3000) | 率 [/周] |
|---|---|---|---|
| A2 の種 11 | (1.81 ± 0.25)e-4 | (3.31 ± 0.47)e-4 | (1.0 ± 0.36)e-7 |
| A2 の種 12 | (1.76 ± 0.22)e-4 | (3.04 ± 0.41)e-4 | (0.86 ± 0.31)e-7 |
| B3 の種 31 | (1.62 ± 0.19)e-4 | (3.20 ± 0.32)e-4 | (1.05 ± 0.25)e-7 |
| B3 の種 32 | (1.96 ± 0.30)e-4 | (3.96 ± 0.50)e-4 | (1.33 ± 0.39)e-7 |
| **4 本 (480 粒子) の重み付き平均** | | | **(1.04 ± 0.16)e-7** |

- 比較 (4.2 GeV・12 MV): 追跡 (10.4 ± 1.6)e-8、線形の式 2.32e-8 (0.22 倍、5σ 低い)、側帯の式 4.46e-8 (0.43 倍、3.7σ 低い)、線形の軌道 + 厳密なスピン回転の MC 5.01e-8 (0.48 倍)。
- **側帯の項で追跡との比は 0.22 → 0.43 に上がったが、追跡の約半分に留まる。MC (線形の軌道) もほぼ同じ 0.48 倍なので、残りの 2 倍は式の導出ではなく、
  線形の軌道で表せない効果 (軌道の非線形: 6 極なし・四極の縦ずれだけのこのリングでは、振幅依存のスピンチューン・x-y-z の非線形結合・ビームの裾の非ガウス分布・RF バケットの非線形) による**
  (どれが主かは確かめていない)。v2 の Polarization は 4.2 GeV・12 MV で脱偏極率を追跡の 0.4 倍に過小評価する。A2 の 0.25 ± 0.06 からは改善 (線形の式のまま)。

### 2.4 大きなリングでの値 (側帯の項と線形の式の比)

| リング | 条件 | 線形 [/周] | 側帯込み [/周] | 比 |
|---|---|---|---|---|
| 大きなリング A (3 GeV、603 m) | 60 台の四極の縦ずれ rms 0.1 mm | 8.79e-14 | 8.84e-14 | 1.006 |
| PF リング (2.5 GeV、187 m) | 同 rms 0.4 mm | 3.45e-9 | 3.50e-9 | 1.014 (P = 0.068: 脱偏極が強い条件) |
| ASTRID2 (0.58 GeV、45.7 m) | 同 rms 1 mm | 1.317e-12 | 1.432e-12 | 1.088 |
| SP8 (大きなリング B、8 GeV、1436 m) | 同 rms 0.1 mm | 未取得 | 未取得 | POL の Emittance[] が 20 分以内に終わらず (スピン行列に NaN の警告)。平面 (ずれなし) の P・τ は A2 の表のとおり |
| 12 セル FODO 3 GeV・1.5 MV | 四極の縦ずれ ±3 mm | 5.52e-11 | 6.26e-11 | 1.13 (A の追跡 3.4e-11、統計 ±22%) |
| 12 セル FODO 2.903 GeV・1.5 MV | 同 | 1.283e-9 | 1.201e-9 | 0.94 (A の追跡 1.05e-9、統計 ±24%) |

- 大きなリング A の rms 0.2・0.4 mm (A2 の §3.3 の表、seed 3、60 台) は **リニア光学が不安定** (1 周写像の x の固有値 1.20 と 0.83 がν = 0.5 に出る。tunes の x が 0.5 になる)。
  A2 の表のこの 2 行 (P = 0.9221、0.9176) は意味のある値ではない。側帯の項は不安定なリングでは使わず線形の式に戻す (`sideband_error`)。0.1 mm は安定で、上の値。
- 計算時間: Sobol 1024 点で、12 セル FODO は +3 秒、大きなリング A は +15 秒 (3563 要素の鎖)。

### 2.5 K64 の `spdepol` との比較

K64 の式は同じ 4.2 GeV・12 MV の例 (`B3_pol_sideband.sad`) で P = 1.265e-4 (gs。pysad k64 も一致)、v2 の線形は P = 0.205、側帯込みは P = 0.119。脱偏極率に直すと
K64 は 4.3e-5 /周 で、独立な値 4.5e-8 の約 1000 倍。K64 の式の位相の広がり ssprd (= 0.056) を 1e-6 まで小さくしても 4.294e-5 のまま (ssprd = 0.056 で 4.300e-5、0.3 倍で 4.294e-5):
**K64 の値はほぼ全部が幅 0 の極限の項** (側帯のガウスの平均の項 `cerfc` はほとんど効かない)。この幅 0 の項は線形の SLIM の率 2.318e-8 の **1852 倍 = 2π × 294.8** で、
A の B55 の指摘 (側帯の項に 2π が余分、1/減衰率の因子) と合う (294.8 は 1/減衰率 (振幅 3.4e-4・5.1e-4・6.0e-4 /周の逆数 2900・2000・1700) の 1/10 程度。厳密な因子は未確認)。
K64 の式の導出に戻らないと直せないので、K64 の式は直さず、独立な導出で置き換えた (`k64sp.py` 相当の実験は §2.5 のこの数値だけ。スクリプトは local/scratch_B3)。

### 2.6 使える範囲と条件 / 使えない範囲

- 使える: RF が有効で放射減衰があり安定な、スピン軌道結合が n の線形近似 + 位相の変調で表せるリング。4.2 GeV の 12 セルリングで 線形の軌道の MC と 10% 以内。FLUC つきの追跡には 0.43 倍 (§2.3) で、脱偏極率は依然過小評価 (P を高く出す)。
  側帯の効果は変調指数 aγσδ/νs ≳ 0.3 で効き始め、0.5〜1 で 1.5〜2.2 倍 (12 セル FODO)。低エネルギー (3 GeV 以下) や大きなリング A・PF では 1〜15% の補正にとどまる。
- 使えない・注意: (1) 結合項 b の u の 2 次の項は入れていない (MC に対し −10%)。(2) 軌道の非線形・振幅依存のスピンチューン・高次の共鳴 (運動論的偏極を含む) は入っていない。
  追跡との比較は §2.3。(3) 結合した光学・νx − νy 近傍など |eig M| が 1 に近い (放射減衰が弱い) と和の項数 K が増え遅い (K は eig M の最大の大きさから決める。最大 3e5 項)。
  (4) RFSW OFF では側帯が定義できず線形の式。(5) 逆向きの磁場のあるリング (MPW) で n が傾く場合 (非平面) の緩和率 Σ r(b·n)² の最適性は A2 から未確認のまま。
- 既存の試験の更新: `test_v2_A_spin.py` の 3 か所 (線形の `depol_rate` は `depol_rate_linear` で比べる、P = 0.8840 → 0.8798)。

## 3. 項目 2: RADPOL の緩和が (1 − PST)/τ の 1.58 倍速い件

### 3.1 K64 の式 (photons.f の `sprot`、`tradk1`・`tradkf1`)

`dst = (PST − st)·cfpd·(p·p0/h²)²·ucde/bt`、スピンに `dst·(bt の向き)` を足す (bt は運動量に垂直な磁場)。`ucde = uc·de`: NOFLUC では `de = cave·anp·uc`、FLUC では `de = dph·uc`。
`uc = cuc·h³/(p0 ρ)` = (3/2) ƛ_C γ²/ρ、`anp = anrad·h·θ`、`cave·anp·uc` = (2/3) r_e γ³ θ/ρ (古典の放射損失)。したがって
`cfpd·uc·de = (5√3/8) r_e ƛ_C γ^5 θ/ρ² = (5√3/8) r_e ƛ_C γ^5 ds/ρ³` で、Sokolov–Ternov の率 1/τ_ST の 1 刻みの分と**厳密に等しい**。
dP = (PST − P)·(1 刻みの率) なので、P=1 からの緩和は (1 − PST)/τ、P=0 からの立ち上がりは PST/τ。`(p p0/h²)² ≈ β⁴` は 1 − 1e-6 以内。

### 3.2 測定

- 3 GeV の 24 セル FODO (τ = 143 s、1 刻みの増分 ~1e-12 ≫ ulp): 追跡 (RADPOL ON、NOFLUC) の P=1 からの 1 周あたりの dP が −1.87320e-10、式 −(1−PST)·sdamp = −1.87318e-10 (9e-6)。
  P=0 からは +2.26965e-9、式 +PST·sdamp = +2.26965e-9 (1.7e-6)。v2・k64 とも (追跡は同じ)。
- ASTRID2 (0.58 GeV、τ = 28815 s、1 周 6336 刻み、1 刻みの (1−PST)·cfpd·ucde = 6.4e-17): 1 周の ucde の和 4.887e-12 に cfpd を掛けると 5.2908e-12 = Emittance[] の sdamp。
  式どおりの増分は ulp (|s|≈1 の下側で 1.1e-16、その半分 5.5e-17) より小さく、倍精度の足し算 s − 6.4e-17 は 1 − 1.1e-16 に丸まる (1.73 倍)。
  P0 = 1、0.5、0.25 から始めた追跡 (1000 周): 1 周あたりの dP と式の比 **1.578、0.934、0.979** (`B3_radpol_roundoff.py`)。P0 = 0.5・0.25 は増分が 3.2・10 ulp で、
  丸めの量子化 (3.2 → 3 ulp = 0.94、10.2 → 10 ulp = 0.98) と一致する。P0 = 1 の 1.58 は 0.58 ulp が 1 ulp に丸まる効果 (1.73) の混ざり。A2 の「1.58 倍」はこれ。
- A2 の記録の「P=1 からの緩和率の絶対値」と「n に垂直な始点からの立ち上がり (式どおり)」が食い違って見えたのも同じ原因 (立ち上がりは |s| ≈ 0 から始まるので丸めが効かない)。
  A2 の `pol4.py` は S_rel = Ssig − b と書いていた (正しくは a − b。比 1.044 になる) が、値の大小の結論は変わらない。

### 3.3 結論と MPW の比

K64 の式に誤りは無く、v2 も追跡を変えない。MPW あり・なしの緩和率の比 (追跡 1.311、v2 1.315、式 1.317) は保たれる (追跡のコードを触っていない)。
低エネルギーのリング (ASTRID2 級) で P≈1 から RADPOL の緩和を測るときは、倍精度の丸めで ST の値が実際とずれる (P0 ≲ 0.5 から始めるか、P=0 からの立ち上がりで測る)。
gs (Fortran の倍精度) も同じ丸めを受ける。

## 4. 項目 3: 光子の記録 (`tphrec`・`tphotonlist`)

### 4.1 実装

- `track6.py`: `_PH` (記録先・幾何・現在の要素の番号 l・`tsetpcvt` 相当の情報)。`_tradk_core` の FLUC の枝で、PHOTONS のとき `tdusrnpl` (tdusrn と同じ乱数の引き方で、光子ごとのエネルギーと位置も返す)
  を使い、光子ごとに `_ph_rec` (tphrec): 放射点 xr = x − r(px − dpx r/2)al、`radangle1` (放射角と偏光。tgauss を 2 回)、`tphchge` (`lattice._tchge_cod` の出口)、
  `tforbitgeo`、`tfgeofrac` (`GeoResult.geo_frac`)。PhotonList の行 = [E(eV)、gx gy gz、ux uy uz、xi1 xi2 xi3、粒子番号、要素番号] (tphotonlist)。
- 対応する要素: BEND (K1=0: `tbend` の刻みごと、K1≠0: `tbendi` の刻みごと + 最後)、QUAD (`tquad`、`tsolqur` の刻み)。K64 と同じ乱数の順序・`pcvt%fr0` の進め方 (tbend.f:367、tbendi.f:141,178、tquad.f:71,108、tsolqu.f:453)。
  MULT・SEXT・SOL 内の要素・ドリフトの放射の光子は位置の情報 (tsetpcvt) を持たないので**記録しない** (乱数だけ K64 と同じに引く)。BEAMBEAM と BEND を同じラティスに置いたときは SADError。
- `ffscmd.py`: BEND があると高速追跡 (numba の TrackRun) をやめ、`track6.track_particles` (Python) で追う。1 粒子 1 周が遅い (24 セルで 1 粒子 1 周 ~1 ms)。

### 4.2 gs との照合 (k64 モード)

`B3_photons_bend.sad`: 3 GeV、2 セル (K1=−0.05 の BEND に DX=1 mm・DY=−0.5 mm・DROTATE=2 mrad、ROTATE=15 度・DX・DY・DROTATE の BEND、普通の BEND、四極)、4 粒子、`SeedRandom[7]`。
光子 358 個 (BEND 6 個で 62、57、67、61、40、70 と QUAD で 1) が、gs と光子数・位置 (GX GY GZ: 全桁)・方向余弦 (10 桁)・偏光 (10 桁)・粒子番号・要素番号まで一致する。
エネルギー列 (E) だけ相対 1.45e-8 (ほぼ一定。多数の光子で 1.40〜1.46e-8、少数の光子で 1.1e-8〜2.8e-8) 違う。E = dp·amass·γ の電子質量の定数 amass の違い (pysad は 0.51099895e6 eV) と推定。確かめていない。

### 4.3 解析式との照合 (v2)

24 セル FODO (3 GeV、ρ = 4.58 m、100 粒子 1 回通過、光子 58 221 個):

| 量 | 式 | 値 | 差 |
|---|---|---|---|
| 光子数 N | (5/(2√3)) α γ θ × 24 磁石 × 100 粒子 = 58 280 | 58 221 | −0.10% (Poisson の 1σ = 0.41%) |
| 臨界エネルギー u_c = (3/2) ħcγ³/ρ | 13 067 eV | (平均 E/u_c から) | — |
| <u>/u_c | 8/(15√3) = 0.3079 | 0.3059 ± 0.0023 | −0.7% (0.9σ) |
| <u²>/u_c² | 11/27 = 0.4074 | 0.3986 | −2.2% (裾の統計) |
| スペクトル (3/(5π)) ∫_y^∞ K_{5/3} | KS 検定 (1% 点 1.63/√N) | 通る | — |
| 位置 | 円弧 (geo_frac) からの距離 | < 2 mm (粒子の軌道のずれ: 分散 × 放射で減った δ) | — |
| 方向 | 接線からの角度 (u > 0.5 u_c の光子) | rms 0.3〜3 /γ | — |

### 4.4 K64 の誤り (v2 で直す): 偏光の xi1

`tphrec` の `c1 = (cost·ppx1 + sint·ppy1)/pp` は、すでに規格化した `ppx1 = ppx/pp`・`ppy1 = ppy/pp` をさらに pp (≈ 曲げ角 数 mrad) で割る。
c1・s1 は偏光の回転角 χ の cos・sin のはずなので |c1| ≤ 1。K64 の xi1 は ±1 を超え 1e3〜1e6 (gs の PhotonList の 8 列目が −3735.37、−1.0015e6 など)。
v2 は pp で割らない: |xi1|、|xi3| ≤ 1、xi1² + xi2² + xi3² ≈ 1 (radangle1 の xi2² + xi3² = 1)。`test_stokes_*`。

### 4.5 使える範囲と条件 / 使えない範囲

- 使える: BEND・QUAD の光子の数・エネルギー・位置・方向・偏光 (RAD・FLUC・PHOTONS ON の TrackParticles)。DX・DY・ROTATE のある BEND も。
- 使えない: MULT・SEXT・SOL 内・ドリフトの光子 (記録されない)。BEAMBEAM と BEND の併用。POL の追跡との併用は未確認。QUAD の光子の位置は gs と未照合 (BEND は全桁)。

## 5. 項目 4: 幾何の fit の DX・DY (tdgeo の qdtrans・tftmatu)

- 状況: SOL 領域 (BOUND≠0 の SOL の組の間) の外の四極の DX・DY は設計の幾何を動かさない (K64 の tdgeo は 0、fit の微分は 0)。領域の中は tsgeo が軌道を追うので動く。
  K64 の tdgeo は、`qdtrans` (要素の出口の軌道の eps = 1e-6 の中心差分) を `tftmatu` (Twiss の写像) で領域の終端まで運び、終端の位置・向きの変化から条件の位置の GX GY GZ CHI1..3 の微分を出す。
- 前提の追加: pysad は領域内の四極のずれ (DX・DY) を、追跡 (`track6n_ext.compile_sol_quad`)・Emittance[] (`solemit._tquade`)・tsgeo (`sgeo._quad_step`) で NotImplementedError にしていた。
  K64 の `tquad` は `tsolrot` で座標を (dx, dy) ずらす。ソレノイドの一様磁場のベクトルポテンシャルのゲージが変わるので、正準運動量もずれる: x → x − dx、y → y − dy で
  px → px + (bz/2)dy、py → py − (bz/2)dx (bz は運動量で規格化、出口は逆)。k_sol_quad に dx, dy を足した (numba のキャッシュを削除)。gs と LINE["GX"…] が全桁一致 (領域内の DX・DY・ROTATE つきの 2 台)。
- 解析微分 (v2): 要素の写像は c_out = F(c − t) + t (t = (DX, −(bz/2)DY, DY, (bz/2)DX))。したがって d c_out/d DX = (I − J) ∂t/∂DX (J = 要素の Jacobian)。
  領域の終端までは要素の写像 (軌道に依存する Jacobian) の積で運ぶ (`match6d.sol_dcod`)。tdgeo の式 (`tdgeo_term`) は K64 の直訳。
- 結果 (`B3_geofit_quad_dx.sad`: QF の DX と QD の DY で M2 の GY = 0.0008、GZ = −0.0035): gs と v2 (解析微分) は 2 回の Newton で収束し、
  QF DX = 2.78848432670e-3 (gs) / 2.78848432669e-3 (v2)、QD DY = 2.66798056318e-3 (両方)。最初の反復の残差は gs 4.86e-8、v2 解析 5.19e-8、v2 中心差分 2.54e-8
  (どの微分も 2 回で収束。K64 と v2 解析の 1 回目の違いは微分の近似の違い。この 1 回目の残差の差は確かめていない)。
- K64 の CHI3 の微分の式は、領域の終端で GCHI3 が 0 に戻される (tsgeo の `tfderotgeo`) ので幾何の変化と合わない (9e-4 を返すが幾何は動かない)。CHI3 を使う fit は K64 でも収束しない。
- 使える範囲: 領域内の四極 (K1≠0・F1=FRINGE=0) の DX・DY。SEXT・OCT・MULT は `_QDX_TYPES` に入れたが、領域内のずれ (追跡の MULT ずれ) は未確認。四極の F1・FRINGE と領域内の BEND・MULT のずれは未対応のまま。
  領域内の四極のずれのある **POL の追跡 (track6 の経路)** は未対応 (NotImplementedError)。k64 モードは SOL 自体のずれ (trots) が未対応のため、ずれた SOL の fit は使えない。

## 6. 項目 5: NAFF (チューンの抽出)

- `pysad/naff.py`: `naff(z, nfreq)` (複素信号の周波数と振幅、Hann 窓、FFT のピーク → Brent 法 (窓つきの内積の最大) → 成分を引いて次、振幅は窓つき最小二乗)、
  `naff_real` (実数の信号)、`tune_xp` (x, px と Twiss から z = x̂ − i p̂)、`tunes_from_coords` (ターン × 6 の座標から νx νy νz)。単位は 1/周回。
- gs の動作 (確認): MAIN の `ON FOURIER; EMIT USE=RING; TRACK …` は `Use EMIT_TANCE or Emittance[] within FFS.` のあと SIGSEGV (tracka → tinip。`work_main_track.md` の結論のまま)。
  FOURIER の出力は gs から取れず、K64 の wjfit の maxord>0 の出力 (作用変数の Fourier のフィット) は比べられない。pysad の MAIN の TRACK は FOURIER を無視する (wtunem の表だけ)。
  wjfit・wfres・wfbin・wiord の直訳は行っていない (使える形 = Python の関数、の要求に NAFF を当てた)。
- 結果: 既知の信号 (単一: 周波数 1e-11、位相 1e-7、2 成分 (間隔 4.2/N): 1e-6、ノイズ SN 14: 1e-5) と、追跡の座標列 (1024 周、x = 2e-7 m √β) のチューンが CALC と
  12 セル FODO 3e-9 (6 極つき: 振幅の 2 乗で動く、1e-5 → 1e-4 で 100 倍)、大きなリング A 6e-10 / 8e-10、PF 1.4e-10 / 1.3e-10、ASTRID2 4.9e-11 / 5.2e-11 (x / y)。
  MAIN の TRACK の wtunem (FFT の放物線補間) は同じ座標列で 4e-5〜2e-4 (`B3_naff_tunes.out.txt`)。
- 使える範囲: 放射を切った (減衰の無い) 追跡の座標列。周波数の間隔が約 2/N 以上の成分。使えない: 放射減衰のある列 (包絡 exp(−λn) で周波数がずれる)、共鳴線の近傍 (間隔 < 1/N)。

## 7. 項目 6: 反転した BEND (P1) と ASTRID2 + MPW の Polarization

A2 の P1 (v2 は緩和率に符号に依らない Σ r_i (b_i·n)² を使う) を、側帯の項を入れた後のコードで確認した (`test_astrid2_mpw_polarization_after_sideband_v2`)。ASTRID2 + MPW (11 周期、逆向きの極):
P = 0.7019、τ_p = 21915 s、独立な式 (∮ρ^-3 の符号つき・絶対値の比) は P = 0.70163、τ = 21886 s で、2e-3 以内 (A2 と同じ)。平面のリング (n が鉛直、ずれなし) では脱偏極率 0 のままで、側帯の項は効かない。
RADINT の K0 の極 (R4、A2) と合わせて、MPW 入りの U0・Jx・τ・P は A2 の範囲で使える (εx は 3.5% の差、A2 の要判断のまま)。確認していない: ずれのある MPW のリング (n が傾く) の P・τ。

## 8. 要判断・残した項目

| 項目 | 見つけたこと | 残した理由 |
|---|---|---|
| 側帯の式と MC の −10% | b(u) の 2 次の項 (変調 × 結合) を入れていない | 導出を 2 次まで広げると M の固有分解の和が 36 倍になる。効果は 10% 以内 |
| 追跡との比較 (§2.3) | 側帯の式は追跡の 0.43 倍 (線形の軌道の MC と同じ 0.48 倍)。残りは軌道の非線形の効果 | 非線形の寄与の分離には、非線形の軌道 + 線形のスピン結合などの切り分けの MC が要る。追跡は 1 粒子 1 周 4.5 ms で、統計誤差 15% (480 粒子) |
| K64 の `spdepol` の誤りの箇所 | 結果が 1000 倍大きい。2π と 1/減衰率の因子 (B55) は未確認 | K64 の式の導出に戻る必要がある。v2 は独立な式で置き換えた |
| 光子のエネルギー列の 1.45e-8 | gs と ほぼ全光子で 1.4e-8 (少数の光子 1.1e-8〜2.8e-8) | amass の定数の差と推定 (gs の値を読んでいない) |
| MULT・SEXT・SOL 内の光子 | 記録しない (K64 は pcvt の値で記録する) | 位置の情報 (tmulta・tmulti の tsetpcvt) の直訳が要る |
| 領域内の四極のずれの POL 追跡 (track6) | NotImplementedError | track6 の `_tquad_sol` に同じずらしを入れれば足りる (未確認) |
| MAIN の TRACK の FOURIER (wjfit の maxord>0) | gs は SIGSEGV で比べられない | 解析的な根拠だけで直訳すると検証できない |
| 大きなリング A の rms 0.2・0.4 mm (A2 §3.3) | リニア光学が不安定 (x の固有値 1.2) | A2 の表の 2 行は意味が無い。0.1 mm のみ使える |

## 9. 走る例 (`example/v2_next/`)

| 例 | 内容 |
|---|---|
| `B3_pol_sideband.sad` (.gs.txt .pysad.txt .pysad_k64.txt) | 4.2 GeV の 12 セル FODO の Polarization。gs・k64: 1.27e-4、v2: 0.119 (側帯込み) |
| `B3_spin_mc.py` | 線形の軌道 + 厳密なスピン回転の MC (§2.2) |
| `B3_radpol_roundoff.py` (.pysad.txt) | ASTRID2 の RADPOL の P0 = 1、0.5、0.25 からの緩和と丸め (§3) |
| `B3_photons_bend.sad` (.gs.txt .pysad.txt .pysad_k64.txt) | 光子 358 個 (gs と k64 が全列一致、v2 は xi1 が ±1 以内) |
| `B3_geofit_quad_dx.sad` (.gs.txt .pysad.txt) | SOL 領域内の四極の DX・DY で幾何の fit (gs と v2 が 1e-12 で一致) |
| `B3_naff_tunes.py` (.out.txt) | NAFF 対 CALC 対 wtunem (4 つのリング) |

`py -3.12 run_sad.py example/v2_next/<名前>.sad` (v2 が既定、`PYSAD_COMPAT=k64` で K64 と同じ値)。gs との並置は `example/tools/run_example.py`。
