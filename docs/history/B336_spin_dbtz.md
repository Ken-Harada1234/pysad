# v2.5.0 記録: Emittance[] (POL) のスピン軌道結合の dbtz の割り算 (B336)

日付: 2026-11-23。試験 `pysad/tests/script_tests/test_v250_spin_dbtz.py` (9 件)。

## 1. 症状

`photons.f` の `tradke` (1537-1552 行) は、スピンの回転ベクトルの磁場の横成分を `btx=xpx/pr**2-blx`、`bty=xpy/pr**2-bly`、`btz=xpzb/pr**2-blz` と作り、座標についての微分 (srot の 4:9 列の元) を `dbtx=dxpx/pr**2-dblx` (`dbtx(6) -= 2*xpx/pr**3`)、`dbty` (同様)、`dbtz=dxpzb/pr-dblz` (`dbtz(6) -= 2*xpzb/pr**3`) とする。z 成分だけ pr で割り、δ の列の補正 (−2 xpzb/pr³) は pr⁻² の微分の形。自分の btz の微分として不整合 (誤植)。pr = 1 + δ (粒子の運動量) なので δ ≠ 0 でだけ効く。

## 2. 確認

gs では確認できない (K64 の Polarization は NaN (B55)、srot の微分列 (4:9) は Emittance[] の SpinPrecessionMatrix (srot(:,1:3)) など他の出力に出ない)。pysad の k64 は誤植をそのまま写している (`spin.py` の元の注記「K64 のまま (pr**2 ではなく pr)」)。

独立な基準: 追跡 (`track6.track_particles` の spins) のスピン回転行列の中心差分から作るスピン軌道結合 D (3×6)。`test_v2_A_spin._spin_orbit_d` と同じ手順。D の解析式は `emit6d.tturne` の `st.srot[:, 3:9]`。3 GeV/c、1 要素 + ドリフト 1 m、軌道 (x, px, y, py) = (8e-3, 1e-3, 6e-3, −1e-3)、δ を変える。|D_解析 − D_追跡| の最大 / |D_追跡| の最大 (v2 の他の修正はそのまま。`spin.DBTZ_PR2` で割り算だけを切り替え):

| 要素 | δ | pr で割る (K64 の式) | pr² で割る (v2) |
|---|---|---|---|
| 矩形の偏向磁石 (15°、L = 1.2 m、K1 = 0) | 0 | 0.80% | 0.80% |
| 同 | 0.05 | 13.2% | 1.20% |
| 同 | 0.10 | 26.9% | 2.24% |
| 同 (K1 = 0.3) | 0 | 0.25% | 0.24% |
| 同 | 0.05 | 6.3% | 0.18% |
| 同 | 0.10 | 12.4% | 0.15% |
| 四極 (L = 0.3、K1 = 0.7) | 0〜0.10 | 0.19% | 0.19% (δ に依らない) |

誤差は δ に比例して増える (13.2% → 26.9%)。z の行の最大誤差は δ = 0.05 の K1 = 0 の偏向磁石で 0.155 (pr) 対 0.014 (pr²) (最大値の相対ではなく絶対)。

## 3. 修正 (v2 のみ)

`spin.tradke_spin` の `dbtz = dxpzb / (pr ** 2 if (_compat.v2() and DBTZ_PR2) else pr) - dblz` (`DBTZ_PR2` は試験が切り替えるための印、既定 True)。k64 は pr。

## 4. 限界・未確認

- 1 要素 + ドリフトで確かめた。リングでの Polarization・脱偏極への影響 (δ のある閉軌道、DP0) は未確認。
- δ = 0.1 の K1 = 0 の偏向磁石の v2 の 2.24% の残り (δ = 0 で 0.8%) は別の項 (端のフリンジ・他の δ 依存。未調査)。
- TrackParticles のスピン追跡 (`track6._sprot`) は別のコードで、この誤植の影響を受けない。
