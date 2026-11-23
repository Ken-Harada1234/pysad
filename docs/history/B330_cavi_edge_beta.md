# v2.5.0 記録: CAVI の端の収束の 1/β (B330)

日付: 2026-11-23。試験 `pysad/tests/script_tests/test_v250_cavi_edge_beta.py` (23 件)。gs の保存した出力 `pysad/tests/script_tests/scripts/v250_cavi_edge_beta__{3gev,2mev,600kev}.{sad,gs.txt}`。確認スクリプト `physics_guide/scripts/cav_02_edge_lorentz.py`。

## 1. 症状

CAVI の端 (入口・出口) の収束と端のエネルギー変化 (`tcav.f` の `tcavfrin`、CALC・Emittance の `tcave.f` の `tcavfrie`) が、ヘルプの H_f = −(eV/L)(sin(ωt − dphi) + …)(x² + y²)/4 を β = 1 のまま使い、粒子の速度 β = p/h を含まない。Δp_x/p0 = x·v/(2 L p0)·(sin(ωt − dphis) + s0)。Maxwell 方程式の場の中のローレンツ力では、端を通る間の運動量変化は ∫F_r dt = ∫F_r ds/(βc) なので Δp_x/p0 = x·v/(2 L p0 β)·(…) (1/β がある)。低エネルギーの空洞 (0.6 MeV/c で β = 0.761) で端の収束が 24% 足りない。

## 2. 確認

### 2.1 gs

L = 1e-4 m、VOLT = 10 V、FREQ = 500 MHz、PHI = 0.3 の CAVI の TrackParticles 1 要素 (FRINGE = 1 と DISFRIN = 1)、x = 1e-5、RFSW ON・RAD OFF。端の収束 (差) と x·v·sin(PHI)/(2 L p0) の比:

| p | β | gs の端の px | 比 (式に 1/β なし) |
|---|---|---|---|
| 3 GeV/c | 1.0000 | 4.925336777684e-11 | 1.000000 |
| 2 MeV/c | 0.9689 | 7.388005165834e-08 | 1.000000 |
| 0.6 MeV/c | 0.7613 | 2.462668381259e-07 | 1.000000 |

gs には 1/β が入っていない。k64 の追跡は gs の px・δ と 12 桁以上一致。

### 2.2 ローレンツ力の数値積分 (独立な基準)

軸の近くの場 (∇·E = 0 と ∇×B = c⁻²∂_tE の展開): qE_z = −(v/L) g(z) sin(ωt − dphis) − (r²/4)(∂_z² − ∂_t²)(同)、E_r = −(r/2)∂_zE_z、B_θ = (r/2c²)∂_tE_z。g = ½(1 + tanh(z/d))、d = 1e-5 m、E0 = 1e5 V/m (L = 2e-4 m、VOLT = 20 V)、500 MHz、電子、DOP853 (rtol 1e-13)。dphis の符号と ΔE の大きさは追跡の x = 0 の ΔE (−v sin(ωt − dphis)) に合わせた (数値積分とのずれ 0.2% 以下)。

入口の px/x (追跡 / 数値積分):

| p | PHI | k64 | v2 |
|---|---|---|---|
| 2 MeV/c (β = 0.9689) | 0.3 | 0.962 | 0.993 |
| 2 MeV/c | 1.2 | 0.968 | 0.999 |
| 0.6 MeV/c (β = 0.7613) | 0.3 | 0.754 | 0.991 |
| 0.6 MeV/c | 1.2 | 0.760 | 0.999 |

k64 の比は β (0.7613 に対し 0.754〜0.760)。端のエネルギー変化の x² の項 (E_z の r² の項まで入れた積分、x = 3e-4): k64 は 0.762 (0.6 MeV/c)・0.970 (2 MeV/c)、v2 は 1.0009・1.0007。δ ≠ 0 (0.6 MeV/c): 数値積分の初期の p_z = p0(1+δ) に対し v2 は δ = −0.05〜+0.15 で 0.991 (一定)、k64 の比は粒子の β (0.7446〜0.8036) に追随する。つまり 1/β は設計の β でなく粒子の β。残りの 0.1〜0.9% は、数値積分が含む本体の Panofsky–Wenzel の項 (O(ωL) = 0.3%) と薄い空洞の近似。出口の端 (FRINGE = 2 と DISFRIN = 1 の差) も同じ (試験は 2% 以内)。

## 3. 修正 (v2 のみ)

- 追跡: `track6._tcavfrin` (python)・`track6n.tcavfrin1` (numba。引数 `b2` を k_cav・k_cav_q・k_cav_in・k_cav_out・cav_frin_all に通す。`track6n._cav_b2()` が `compat.v2()`、`track6n_ext` の呼び出しも同じ)。
- CALC・Emittance: `emit6d._tcavfrie_map` (`solemit._cav_fringe`・`emit6d._tcave` から使う)、`trptacc._cavfrie_piece` (v2 は同じ写像の数値微分。K64 の解析ヤコビアンは 1/β の無い式なので使わない)。
- 式: ib = h1/p1 (1/β、粒子の β) として、kick は Δ(p_x/p0) = x·vf·(sin(ωt − dphis) + s0)·ib、エネルギー変化は dh に ib をかける。H_f ∝ 1/β(E) を一貫させるため、時間のずれ Δt = −∂H_f/∂E = −(x² + y²)·κ0/(2 p1³) (κ0 = vf·p0·(sin(ωt − dphis) + s0)、∂(1/β)/∂h = −1/p³) を入れ、z = −(t + Δt)·p2/h2 にする。時間のずれが無いと、ヤコビアンのシンプレクティック条件の誤差が 0.6 MeV/c・L = 1 mm・VOLT = 1 kV・x = 3 mm・δ = 0.03 で 7e-4 に悪化する (入れると 7.8e-5。k64 は 6e-5 で、主因は端以外)。
- 3 GeV/c (1/β − 1 = 1.5e-8) では差が相対 1e-7 以下。DISFRIN = 1 と L = 0 は変わらない。

## 4. 限界・未確認

- 端 1 つの 1 次の収束とエネルギー変化の x² の項だけを数値積分と比べた。長い空洞 (位相の進み ωL/β が大きい) の本体との組合せは見ていない。K64 の CAVI は、基準粒子の位相が入口・出口・中央で同じ (進行波が基準粒子と同期する) 扱いで、定在波の位相の進みは入れない。
- 時間のずれの項 (x² の 2 次) は、ローレンツ力の積分では確かめられない (r³ までの場が要る)。シンプレクティック条件だけで決めた。
- ヤコビアン (CALC・Emittance) は TRPT の解析式を使わず数値微分。Emittance[] は gs と比べていない (gs は k64 と同じ 1/β の無い式なので、v2 は gs と違う値になる。k64 は従来どおり)。
- 軸の近くの場 (r の 1 次) の積分。TCAVI の端、MULT に RF を持たせた tmultiacc の端 (同じ `tcavfrin` を呼ぶ) は、同じ修正が入るが個別には確かめていない。
